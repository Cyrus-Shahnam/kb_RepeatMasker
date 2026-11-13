#!/usr/bin/env python3
"""
Dry-run the whole Runner OUTSIDE Docker with a fake RepeatMasker + famdb.py and
fake KBase clients (lessons §8).  Exercises param handling, command building,
contig-ID renaming, PATH injection, parsing, warnings, zips and the HTML report.

  python3 tools/dry_run.py            # from the module root; needs only python3
"""
import json
import os
import shutil
import sys
import tempfile
import traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRATCH = tempfile.mkdtemp(prefix="kb_rm_dryrun_")
FAKE_BIN = os.path.join(ROOT, "tools", "fake_bin")
for f in os.listdir(FAKE_BIN):
    os.chmod(os.path.join(FAKE_BIN, f), 0o755)
os.environ["REPEATMASKER_ENV_BIN"] = FAKE_BIN
os.environ["RM_FAMDB_DIR"] = os.path.join(SCRATCH, "data", "famdb")
os.makedirs(os.environ["RM_FAMDB_DIR"])
open(os.path.join(os.environ["RM_FAMDB_DIR"], "dfam40.0.h5"), "w").write("fake")
sys.path[:0] = [os.path.join(ROOT, "lib"), os.path.join(ROOT, "test"),
                os.path.join(ROOT, "tools")]

from kb_repeatmasker.repeatmasker_utils import RepeatMaskerRunner, validate_params  # noqa
import make_fixtures  # noqa: E402
import score_masking  # noqa: E402
from plant_repeats import read_fasta, write_fasta, random_seq  # noqa: E402


class FakeAU:
    def __init__(self, files):
        self.files, self.saved, self.n = dict(files), {}, 0

    def get_assembly_as_fasta(self, p):
        path = self.files[p["ref"]]
        dst = os.path.join(SCRATCH, "dl_%d_%s" % (len(os.listdir(SCRATCH)),
                                                  os.path.basename(path)))
        shutil.copy(path, dst)
        return {"path": dst, "assembly_name": os.path.basename(path).rsplit(".", 1)[0]}

    def save_assembly_from_fasta(self, p):
        assert os.path.dirname(os.path.abspath(p["file"]["path"])) == SCRATCH, \
            "saved FASTA must be in scratch (AssemblyUtil only mounts scratch)"
        self.n += 1
        ref = "99/%d/1" % self.n
        self.saved[ref] = (p["assembly_name"], p["file"]["path"])
        self.files[ref] = p["file"]["path"]
        return ref


class FakeReport:
    def __init__(self):
        self.last = None

    def create_extended_report(self, p):
        for fl in p["file_links"]:
            assert os.path.exists(fl["path"]), fl
        assert os.path.exists(os.path.join(p["html_links"][0]["path"], "index.html"))
        self.last = p
        return {"name": p["report_object_name"], "ref": "99/100/1"}


def runner(au, rep):
    return RepeatMaskerRunner({"scratch": SCRATCH}, None, {"token": "x"},
                              clients={"au": au, "report": rep})


RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append((name, bool(cond)))
    print("%s  %s %s" % ("PASS" if cond else "FAIL", name, detail))


def main():
    fx, man = make_fixtures.build(os.path.join(SCRATCH, "fixtures"))
    au = FakeAU({"1/1/1": fx["planted"], "1/2/1": fx["control"], "1/3/1": fx["library"]})
    rep = FakeReport()

    # 1. custom library, two inputs, soft masking
    r = runner(au, rep)
    out = r.run({"workspace_name": "ws", "assembly_refs": ["1/1/1", "1/2/1", "1/1/1"],
                 "custom_library_ref": "1/3/1", "species": "", "div": "",
                 "cutoff": "225.0", "norna": "0", "output_suffix": "_masked"})
    check("custom: two outputs (duplicate input dropped)",
          len(out["masked_assembly_refs"]) == 2, out["masked_assembly_refs"])
    planted_masked = au.saved[out["masked_assembly_refs"][0]][1]
    sc = score_masking.score(man, read_fasta(planted_masked))
    check("custom: TE copies detected", sc["te_detected"] >= 12,
          "%d/%d recall %.3f" % (sc["te_detected"], sc["te_copies"], sc["te_bp_recall"]))
    check("custom: microsatellites masked",
          all(s["masked_fraction"] > 0.9 for s in sc["simple"]))
    ids_in = [i for i, _ in read_fasta(fx["planted"])]
    ids_out = [i for i, _ in read_fasta(planted_masked)]
    check("IDs round-trip incl. >50-char ID", ids_in == ids_out, ids_out[0][:30])
    rows = r.results[0]["rows"]
    cats = set(r.results[0]["categories"])
    check(".out classes parsed (LINE LTR DNA SINE Simple_repeat)",
          {"LINE", "LTR", "DNA", "SINE", "Simple_repeat"} <= cats, sorted(cats))
    check(".out repeat names are family names, not strands",
          not any(x["repeat"] in ("+", "C") for x in rows))
    check("simple + low-complexity bp counted", r.results[0]["simple_lowcomp_bp"] >= 216)
    check(".out rows carry original IDs",
          any(x["contig"] == make_fixtures.LONG_ID for x in rows))
    check("masked_bp == lowercase count", r.results[0]["masked_bp"] ==
          sum(sum(c.islower() for c in s) for _, s in read_fasta(planted_masked)))
    check("cmd: input last + absolute", r.commands[-1].split()[-1].startswith("/"))
    check("cmd: -lib passed, no -species", "-lib" in r.commands[-1]
          and "-species" not in r.commands[-1])
    check("report: summary TSV + 2 zips", len(rep.last["file_links"]) == 3)
    html = open(os.path.join(rep.last["html_links"][0]["path"], "index.html")).read()
    check("report: HTML has commands + citation",
          "Commands executed" in html and "Dfam" in html)
    check("report: duplicate warning shown", "Duplicate input" in html)
    zips = [fl["path"] for fl in rep.last["file_links"] if fl["path"].endswith(".zip")]
    import zipfile
    names = zipfile.ZipFile(zips[0]).namelist()
    check("zip: .out .gff .tbl masked.fa id map",
          all(any(n.endswith(s) for n in names) for s in
              (".out", ".gff", ".tbl", ".masked.fa", "contig_id_map.tsv")), len(names))
    gff = [n for n in names if n.endswith(".gff")][0]
    gtxt = zipfile.ZipFile(zips[0]).read(gff).decode()
    check("gff uses original IDs", make_fixtures.LONG_ID in gtxt and "rmq0" not in gtxt)

    # 2. species (Dfam) mode, hard N, interspersed only
    os.environ["FAKE_SPECIES_LIB"] = fx["library"]
    os.environ["FAKE_FAMDB_MISSING"] = "1"
    r = runner(au, rep)
    out = r.run({"workspace_name": "ws", "assembly_refs": ["1/1/1"],
                 "species": "Drosophila melanogaster", "mask_type": "hard_N",
                 "mask_scope": "interspersed_only", "cutoff": "300",
                 "uncurated": 1, "speed": "quick"})
    m2 = au.saved[out["masked_assembly_refs"][0]][1]
    sc2 = score_masking.score(man, read_fasta(m2))
    check("dfam: TE copies masked with N", sc2["te_detected"] >= 12)
    check("dfam: -nolow leaves microsatellites", all(s["masked_fraction"] < 0.1
                                                     for s in sc2["simple"]))
    cmd = r.commands[-1]
    check("dfam: -species/-uncurated/-nolow/-q, no -cutoff",
          "-species" in cmd and "-uncurated" in cmd and "-nolow" in cmd
          and " -q " in cmd and "-cutoff" not in cmd, cmd)
    html = open(os.path.join(rep.last["html_links"][0]["path"], "index.html")).read()
    check("dfam: cutoff-dropped warning", "only applies to custom-library" in html)
    check("dfam: missing-partition warning", "dfam40.uncurated.consensus.9.h5" in html)
    check("dfam: famdb validation command recorded",
          any("famdb.py" in c and "lineage" in c for c in r.commands))

    # 3. combined mode
    r = runner(au, rep)
    out = r.run({"workspace_name": "ws", "assembly_refs": ["1/2/1"],
                 "species": "Drosophila melanogaster", "custom_library_ref": "1/3/1",
                 "save_masked_assembly": 0})
    check("combined: no assembly saved when not requested",
          out["masked_assembly_refs"] == [])
    check("combined: -lib with combined library",
          "combined_library.fa" in r.commands[-1])

    # 4. unknown species -> readable error before RepeatMasker runs
    try:
        runner(au, rep).run({"workspace_name": "ws", "assembly_refs": ["1/1/1"],
                             "species": "Notarealus taxonus"})
        check("bad species raises", False)
    except ValueError as e:
        check("bad species raises readable error", "not found" in str(e), str(e)[:70])

    # 5. nothing found -> no .masked file path
    lib2 = os.path.join(SCRATCH, "unrelated_lib.fa")
    import random
    write_fasta(lib2, [("junk#Unknown", random_seq(random.Random(1), 800))])
    au.files["1/4/1"] = lib2
    r = runner(au, rep)
    out = r.run({"workspace_name": "ws", "assembly_refs": ["1/2/1"],
                 "custom_library_ref": "1/4/1", "mask_scope": "interspersed_only"})
    m5 = au.saved[out["masked_assembly_refs"][0]][1]
    check("no repeats: output == input sequence",
          read_fasta(m5) == read_fasta(fx["control"]))
    check("no repeats: warning", any("No repeats" in w for w in r.results[0]["warnings"]))

    # 6. fast param validation
    for bad, frag in (({"assembly_refs": ["1/1/1"], "species": "x y"}, "workspace"),
                      ({"workspace_name": "w", "assembly_refs": [], "species": "Fungi"},
                       "assembly_refs"),
                      ({"workspace_name": "w", "assembly_refs": ["a"]}, "Neither"),
                      ({"workspace_name": "w", "assembly_refs": ["a"],
                        "species": "dros; rm -rf /"}, "unsupported"),
                      ({"workspace_name": "w", "assembly_refs": ["a"], "species": "Fungi",
                        "div": "150"}, "div"),
                      ({"workspace_name": "w", "assembly_refs": ["a"], "species": "Fungi",
                        "mask_type": "lower"}, "mask_type"),
                      ({"workspace_name": "w", "assembly_refs": ["a"], "species": "Fungi",
                        "output_suffix": "a b"}, "output_suffix")):
        try:
            validate_params(bad)
            check("validation rejects (%s)" % frag, False)
        except ValueError as e:
            check("validation rejects (%s)" % frag, frag in str(e))
    m, _w = validate_params({"workspace_name": "w", "assembly_refs": "1/1/1",
                             "species": " Fungi ", "div": "10.0", "cutoff": "",
                             "norna": "1", "mask_type": "", "speed": None})
    check("coercion: blanks->defaults, strings->numbers",
          m["div"] == 10.0 and m["cutoff"] == 225 and m["norna"] == 1 and
          m["mask_type"] == "soft" and m["speed"] == "default" and
          m["species"] == "Fungi" and m["assembly_refs"] == ["1/1/1"])

    n_fail = sum(1 for _, ok in RESULTS if not ok)
    print("\n%d checks, %d failed. Scratch: %s" % (len(RESULTS), n_fail, SCRATCH))
    print("Report to open in a browser: %s"
          % os.path.join(rep.last["html_links"][0]["path"], "index.html"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(2)
