# -*- coding: utf-8 -*-
"""
kb_repeatmasker tests.

Fast tier (seconds, no tool run):  test_fast_*
Known-truth tier (real RepeatMasker, run once per class and cached): test_kt*
  kt0-kt5  synthetic TEs planted into random background, custom library
  kt6      hard-N + interspersed-only run
  kt7-kt8  Dfam species path (needs /data/famdb, else SKIPPED - read the summary!)
"""
import json
import os
import random
import shutil
import sys
import time
import unittest
from configparser import ConfigParser

import requests

from kb_repeatmasker.kb_repeatmaskerImpl import kb_repeatmasker
from kb_repeatmasker.kb_repeatmaskerServer import MethodContext
from kb_repeatmasker import rm_parsers as P
from kb_repeatmasker import repeatmasker_utils as U
from installed_clients.AssemblyUtilClient import AssemblyUtil

try:
    from installed_clients.authclient import KBaseAuth as _KBaseAuth
except ImportError:
    from kb_repeatmasker.authclient import KBaseAuth as _KBaseAuth

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "tools")]
import make_fixtures  # noqa: E402
import score_masking  # noqa: E402
from plant_repeats import plant, random_seq, read_fasta, write_fasta  # noqa: E402

DMEL = "Drosophila melanogaster"


class kb_repeatmaskerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get("KB_AUTH_TOKEN", None)
        config = ConfigParser()
        config.read(os.environ.get("KB_DEPLOYMENT_CONFIG", None))
        cls.cfg = dict(config.items("kb_repeatmasker"))
        user_id = _KBaseAuth(cls.cfg["auth-service-url"]).get_user(token)
        cls.token = token
        cls.ctx = MethodContext(None)
        cls.ctx.update({"token": token, "user_id": user_id,
                        "provenance": [{"service": "kb_repeatmasker",
                                        "method": "please_never_use_it_in_production",
                                        "method_params": []}],
                        "authenticated": 1})
        cls.wsURL = cls.cfg["workspace-url"]
        cls.serviceImpl = kb_repeatmasker(cls.cfg)
        cls.scratch = cls.cfg["scratch"]
        cls.callback_url = os.environ["SDK_CALLBACK_URL"]
        cls.au = AssemblyUtil(cls.callback_url, token=token)
        cls.wsName = "test_kb_repeatmasker_%d" % int(time.time() * 1000)
        cls._ws_call("create_workspace", {"workspace": cls.wsName})
        cls.kt = {}

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "wsName"):
            try:
                cls._ws_call("delete_workspace", {"workspace": cls.wsName})
                print("Test workspace was deleted")
            except Exception as e:  # never mask a test failure with cleanup
                print("could not delete workspace: %s" % e)

    @classmethod
    def _ws_call(cls, method, param):
        r = requests.post(cls.wsURL, headers={"Authorization": cls.token},
                          json={"version": "1.1", "id": "1",
                                "method": "Workspace." + method, "params": [param]})
        r.raise_for_status()
        body = r.json()
        if "error" in body:
            raise RuntimeError(body["error"])
        res = body.get("result")  # JSON-RPC 1.1: a list holding the return value
        return res[0] if isinstance(res, list) and res else res

    @classmethod
    def stage(cls, path, name):
        """AssemblyUtil mounts scratch, not test/: copy the fixture there first."""
        dst = os.path.join(cls.scratch, "%s_%d.fa" % (name, int(time.time() * 1e6)))
        shutil.copy(path, dst)
        return cls.au.save_assembly_from_fasta({"file": {"path": dst},
                                                "workspace_name": cls.wsName,
                                                "assembly_name": name})

    def run_app(self, **params):
        params.setdefault("workspace_name", self.wsName)
        return self.serviceImpl.run_repeatmasker(self.ctx, params)[0]

    # ================================================================ fast tier
    def test_fast_validation_errors(self):
        cases = [({"assembly_refs": ["1/2/3"], "species": "Fungi"}, "workspace_name"),
                 ({"workspace_name": "w", "assembly_refs": [""], "species": "Fungi"},
                  "assembly_refs"),
                 ({"workspace_name": "w", "assembly_refs": ["1/2/3"]}, "Neither"),
                 ({"workspace_name": "w", "assembly_refs": ["1/2/3"],
                   "species": "x; rm -rf /"}, "unsupported characters"),
                 ({"workspace_name": "w", "assembly_refs": ["1/2/3"], "species": "Fungi",
                   "div": "0"}, "div"),
                 ({"workspace_name": "w", "assembly_refs": ["1/2/3"], "species": "Fungi",
                   "cutoff": "abc"}, "cutoff"),
                 ({"workspace_name": "w", "assembly_refs": ["1/2/3"], "species": "Fungi",
                   "speed": "warp"}, "speed"),
                 ({"workspace_name": "w", "assembly_refs": ["1/2/3"], "species": "Fungi",
                   "output_suffix": "has space"}, "output_suffix")]
        for params, frag in cases:
            with self.assertRaises(ValueError) as cm:
                U.validate_params(params)
            self.assertIn(frag, str(cm.exception))

    def test_fast_narrative_blank_and_string_params(self):
        m, w = U.validate_params({
            "workspace_name": "w", "assembly_refs": "1/2/3", "species": "Fungi",
            "custom_library_ref": "", "div": "15.0", "cutoff": "250.0", "norna": "1",
            "uncurated": "0", "mask_type": "", "mask_scope": None, "speed": "",
            "alignments": "", "save_masked_assembly": "1", "output_suffix": ""})
        self.assertEqual(m["assembly_refs"], ["1/2/3"])
        self.assertEqual(m["div"], 15.0)
        self.assertEqual(m["cutoff"], 225)  # dropped: Dfam mode
        self.assertTrue(any("only applies to custom-library" in x for x in w))
        self.assertEqual((m["norna"], m["uncurated"], m["mask_type"], m["speed"]),
                         (1, 0, "soft", "default"))
        self.assertEqual(m["output_suffix"], "_masked")
        self.assertEqual(m["library_mode"], "dfam")

    def test_fast_command_order(self):
        m, _ = U.validate_params({"workspace_name": "w", "assembly_refs": ["a"],
                                  "custom_library_ref": "b", "cutoff": 300,
                                  "mask_type": "hard_X", "speed": "slow"})
        cmd = U.build_rm_command("/x/RepeatMasker", "rel/q.fa", "out", m, 2, "lib.fa")
        self.assertEqual(cmd[-1], os.path.abspath("rel/q.fa"))  # input LAST, absolute
        self.assertTrue(os.path.isabs(cmd[cmd.index("-dir") + 1]))
        self.assertTrue(os.path.isabs(cmd[cmd.index("-lib") + 1]))
        for flag in ("-x", "-s", "-gff", "-excln"):
            self.assertIn(flag, cmd)
        self.assertEqual(cmd[cmd.index("-cutoff") + 1], "300")
        self.assertNotIn("-species", cmd)
        self.assertNotIn("-xsmall", cmd)

    def test_fast_parsers(self):
        d = os.path.join(self.scratch, "parser_fixture")
        os.makedirs(d, exist_ok=True)
        out = os.path.join(d, "q.out")
        with open(out, "w") as fh:
            fh.write("   SW   perc perc perc  query    position in query   matching\n"
                     "score   div. del. ins.  sequence  begin end (left)  repeat\n\n"
                     "  463  12.1  0.9  1.2  rmq000001     101   400 (69600) + "
                     "kbTE_SINE1 SINE/tRNA  1  300  (20)   1\n"
                     "  255   0.0  0.0  0.0  rmq000002    5001  5120 (  44880) + "
                     "(CA)n      Simple_repeat  1  120  (0)   2 *\n"
                     " 2210   5.0  1.0  1.0  rmq000002    9000 10499 (39501) C "
                     "kbTE_LTR1  LTR/Gypsy (0) 1500  1   3\n")
        idmap = {"rmq000001": make_fixtures.LONG_ID, "rmq000002": "ctg2"}
        rows, none = P.parse_out(out, idmap)
        self.assertFalse(none)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["contig"], make_fixtures.LONG_ID)
        self.assertEqual(rows[2]["strand"], "-")
        self.assertTrue(rows[1]["overlapped"])
        s = P.summarize_rows(rows)
        self.assertEqual(s["interspersed_bp"], 300 + 1500)
        self.assertEqual(s["simple_lowcomp_bp"], 120)
        tr = os.path.join(d, "t.out")
        P.translate_out_file(out, tr, idmap)
        txt = open(tr).read()
        self.assertIn(make_fixtures.LONG_ID, txt)
        self.assertNotIn("rmq000001", txt)
        none_out = os.path.join(d, "none.out")
        open(none_out, "w").write("There were no repetitive sequences detected in q.fa\n")
        self.assertEqual(P.parse_out(none_out), ([], True))
        self.assertEqual(P.parse_famdb_check(
            "  Uncurated HMMs  partition 50 [Eutheria]:  MISSING  [dfam40.uncurated.hmm"
            ".50.h5]"), ["dfam40.uncurated.hmm.50.h5"])
        self.assertIsNone(P.parse_lineage_totals("No species found for search term 'x'"))

    def test_fast_fixture_sanity(self):
        paths, man = make_fixtures.build(os.path.join(self.scratch, "fx_sanity"))
        planted = read_fasta(paths["planted"])
        control = read_fasta(paths["control"])
        self.assertEqual([c["id"] for c in man["contigs"]], [i for i, _ in planted])
        self.assertEqual([c["length"] for c in man["contigs"]],
                         [len(s) for _, s in planted])
        self.assertEqual([len(s) for _, s in planted], [len(s) for _, s in control])
        self.assertGreater(len(make_fixtures.LONG_ID), 50)
        self.assertEqual(sum(p["kind"] == "te" for p in man["plants"]), 16)
        # planting is the ONLY difference between planted and control
        diff = sum(a != b for (_, x), (_, y) in zip(planted, control)
                   for a, b in zip(x, y))
        self.assertLessEqual(diff, sum(p["length"] for p in man["plants"]))

    # =========================================================== known truth
    @classmethod
    def lib_run(cls):
        """Planted + control in ONE job with the synthetic library (soft mask)."""
        if "lib" not in cls.kt:
            paths, man = make_fixtures.build(os.path.join(cls.scratch, "fx_kt"))
            refs = {k: cls.stage(paths[k], "kt_" + k)
                    for k in ("planted", "control", "library")}
            impl_out = cls.serviceImpl.run_repeatmasker(cls.ctx, {
                "workspace_name": cls.wsName,
                "assembly_refs": [refs["planted"], refs["control"]],
                "custom_library_ref": refs["library"], "species": "",
                "mask_type": "soft", "output_suffix": "_masked"})[0]
            masked = []
            for ref in impl_out["masked_assembly_refs"]:
                got = cls.au.get_assembly_as_fasta({"ref": ref})
                masked.append(read_fasta(got["path"]))
            cls.kt["lib"] = {"paths": paths, "man": man, "refs": refs,
                             "out": impl_out, "masked": masked}
            print("KT lib run:", json.dumps(impl_out["repeatmasker_results"], indent=1))
        return cls.kt["lib"]

    def test_kt0_outputs_and_report(self):
        k = self.lib_run()
        out = k["out"]
        self.assertEqual(len(out["masked_assembly_refs"]), 2)
        self.assertEqual(len(out["repeatmasker_results"]), 2)
        self.assertTrue(out["report_ref"])
        rep = self._ws_call("get_objects2", {"objects": [{"ref": out["report_ref"]}]})
        data = rep["data"][0]["data"]
        self.assertEqual(len(data["objects_created"]), 2)
        self.assertEqual(len(data["file_links"]), 3)  # summary TSV + 2 zips

    def test_kt1_planted_repeats_are_masked(self):
        k = self.lib_run()
        sc = score_masking.score(k["man"], k["masked"][0])
        print("KT1 planted:", sc["te_detected"], "/", sc["te_copies"], "recall",
              sc["te_bp_recall"], "by divergence", sc["detected_by_divergence"])
        # Tightened after the first real run (2026-11-13: 16/16, recall 0.9994).
        self.assertGreaterEqual(sc["te_detected"], 15)       # of 16 copies
        self.assertGreaterEqual(sc["te_bp_recall"], 0.95)
        for s in sc["simple"]:
            self.assertGreater(s["masked_fraction"], 0.5, s["family"])

    def test_kt2_control_has_no_te_hits(self):
        k = self.lib_run()
        ctrl = k["out"]["repeatmasker_results"][1]
        te_bp = sum(v for c, v in ctrl["categories_bp"].items()
                    if c not in P.NON_INTERSPERSED)
        self.assertLessEqual(te_bp, 200, ctrl["categories_bp"])
        total = sum(len(s) for _, s in k["masked"][1])
        self.assertLess(ctrl["masked_bp"], 0.005 * total)

    def test_kt3_background_not_masked_in_planted(self):
        k = self.lib_run()
        sc = score_masking.score(k["man"], k["masked"][0])
        self.assertLess(sc["background_masked_frac"], 0.005, sc["background_masked_bp"])

    def test_kt4_contig_ids_and_lengths_round_trip(self):
        k = self.lib_run()
        for src, masked in ((k["paths"]["planted"], k["masked"][0]),
                            (k["paths"]["control"], k["masked"][1])):
            orig = read_fasta(src)
            self.assertEqual([i for i, _ in orig], [i for i, _ in masked])
            self.assertEqual([len(s) for _, s in orig], [len(s) for _, s in masked])
        self.assertEqual(k["masked"][0][0][0], make_fixtures.LONG_ID)

    def test_kt5_soft_mask_survives_assembly_save(self):
        k = self.lib_run()
        lower = sum(sum(c.islower() for c in s) for _, s in k["masked"][0])
        reported = k["out"]["repeatmasker_results"][0]["masked_bp"]
        self.assertGreater(reported, 0)
        self.assertEqual(lower, reported,
                         "AssemblyUtil changed case: soft-masked lowercase did not "
                         "survive save/download (%d vs %d). Default to hard masking."
                         % (lower, reported))

    def test_kt6_hard_n_interspersed_only(self):
        k = self.lib_run()
        out = self.run_app(assembly_refs=[k["refs"]["planted"]],
                           custom_library_ref=k["refs"]["library"],
                           mask_type="hard_N", mask_scope="interspersed_only",
                           output_suffix="_hardN")
        got = self.au.get_assembly_as_fasta({"ref": out["masked_assembly_refs"][0]})
        masked = read_fasta(got["path"])
        sc = score_masking.score(k["man"], masked)
        self.assertGreaterEqual(sc["te_detected"], 15)
        for s in sc["simple"]:  # -nolow: microsatellites left alone
            self.assertLess(s["masked_fraction"], 0.2, s["family"])
        n_count = sum(s.upper().count("N") for _, s in masked)
        self.assertEqual(n_count, out["repeatmasker_results"][0]["masked_bp"])

    # ------------------------------------------------- Dfam (reference data)
    def _need_famdb(self):
        if not U.RepeatMaskerRunner({"scratch": self.scratch}).famdb_ready():
            self.skipTest("No Dfam reference data at %s (did init run?)"
                          % U.FAMDB_DATA_DIR)

    def test_kt7_dfam_species_known_truth(self):
        self._need_famdb()
        work = os.path.join(self.scratch, "fx_dfam")
        os.makedirs(work, exist_ok=True)
        runner = U.RepeatMaskerRunner({"scratch": self.scratch})
        lib = os.path.join(work, "dmel_curated.fa")
        runner.famdb(["families", "-f", "fasta_name", "--include-class-in-name", "-ad",
                      "--curated", DMEL], work, os.path.join(work, "famdb.log"),
                     stdout_path=lib)
        fams = [(n, s) for n, s in read_fasta(lib)
                if 800 <= len(s) <= 5000 and any("#%s/" % c in n
                                                 for c in ("LTR", "LINE", "DNA"))]
        if len(fams) < 6:
            self.skipTest("only %d suitable Dfam families for %s" % (len(fams), DMEL))
        random.Random("dfam-kt").shuffle(fams)
        fams = fams[:6]
        bg = [("dmel_bg1", random_seq(random.Random("dfam-bg1"), 60000)),
              ("dmel_bg2", random_seq(random.Random("dfam-bg2"), 40000))]
        planted, man = plant(bg, fams, 6, (0.03,), "dfam-kt")
        fa = os.path.join(work, "dmel_planted.fa")
        write_fasta(fa, planted)
        ref = self.stage(fa, "kt_dmel_planted")
        out = self.run_app(assembly_refs=[ref], species=DMEL, output_suffix="_dfam")
        got = self.au.get_assembly_as_fasta({"ref": out["masked_assembly_refs"][0]})
        sc = score_masking.score(man, read_fasta(got["path"]))
        print("KT7 Dfam %s: %d/%d detected, recall %.3f" % (
            DMEL, sc["te_detected"], sc["te_copies"], sc["te_bp_recall"]))
        self.assertGreaterEqual(sc["te_detected"], 5)        # 2026-11-13: 6/6
        self.assertGreaterEqual(sc["te_bp_recall"], 0.90)    # 2026-11-13: 1.000
        self.assertLess(sc["background_masked_frac"], 0.01)

    def test_kt8_unknown_species_is_readable_error(self):
        self._need_famdb()
        k = self.lib_run()
        with self.assertRaises(ValueError) as cm:
            self.run_app(assembly_refs=[k["refs"]["control"]],
                         species="Notarealus taxonus")
        self.assertIn("not found", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
