"""
RepeatMaskerRunner: all tool logic for the kb_repeatmasker app.

The SDK-managed Impl only parses params and delegates here (lessons §1), so
`make` can never clobber this code.
"""
import collections
import datetime
import glob
import html
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import uuid

from kb_repeatmasker import rm_parsers as P

ENV_BIN = os.environ.get("REPEATMASKER_ENV_BIN", "/opt/conda/envs/repeatmasker/bin")
FAMDB_DATA_DIR = os.environ.get("RM_FAMDB_DIR", "/data/famdb")
REFDATA_ROOT = os.path.dirname(FAMDB_DATA_DIR)
EXPECTED_RM_VERSION = "4.2.4"
EXPECTED_DFAM_VERSION = "4.0"
DEFAULT_CUTOFF = 225

MASK_TYPES = ("soft", "hard_N", "hard_X")
MASK_SCOPES = ("all", "interspersed_only", "simple_only")
SPEEDS = ("default", "quick", "rush", "slow")
SPECIES_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _.'()\-]{1,99}$")
SUFFIX_RE = re.compile(r"^[A-Za-z0-9|._\-]{0,40}$")

DEFAULTS = {
    "species": "",
    "custom_library_ref": "",
    "mask_type": "soft",
    "mask_scope": "all",
    "speed": "default",
    "norna": 0,
    "uncurated": 0,
    "div": None,
    "cutoff": DEFAULT_CUTOFF,
    "alignments": 0,
    "save_masked_assembly": 1,
    "output_suffix": "_masked",
}

CITATIONS = [
    "Smit AFA, Hubley R, Green P. RepeatMasker Open-4.0. 2013-2015. "
    "http://www.repeatmasker.org (Open Software License v2.1)",
    "Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF. The Dfam community resource "
    "of transposable element families, sequence models, and genome annotations. "
    "Mobile DNA 12:2 (2021). doi:10.1186/s13100-020-00230-y (Dfam data: CC0 1.0)",
    "Benson G. Tandem repeats finder: a program to analyze DNA sequences. "
    "Nucleic Acids Res 27:573-580 (1999).",
]


def _blank(v):
    return v is None or (isinstance(v, str) and v.strip() == "")


def _as_int(v, name):
    try:
        return int(float(str(v).strip()))
    except (TypeError, ValueError):
        raise ValueError("Parameter '%s' must be a whole number, got %r" % (name, v))


def _as_float(v, name):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        raise ValueError("Parameter '%s' must be a number, got %r" % (name, v))


def _as_bool(v):
    if isinstance(v, str):
        return v.strip().lower() in ("1", "true", "yes", "on", "1.0")
    return bool(v)


def safe_name(s, maxlen=60):
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_.") or "assembly"
    return s[:maxlen]


def validate_params(params):
    """
    Validate and normalise app parameters before any I/O.
    Returns (merged_params, warnings).  Raises ValueError with a readable message.
    """
    if not isinstance(params, dict):
        raise ValueError("params must be a dict")
    warnings = []
    ws = params.get("workspace_name")
    if _blank(ws):
        raise ValueError("Parameter 'workspace_name' is required")

    refs = params.get("assembly_refs")
    if isinstance(refs, str):
        refs = [refs]
    refs = [r.strip() for r in (refs or []) if not _blank(r)]
    if not refs:
        raise ValueError("Parameter 'assembly_refs' needs at least one Assembly")
    deduped = list(collections.OrderedDict.fromkeys(refs))
    if len(deduped) < len(refs):
        warnings.append("Duplicate input assemblies were removed (%d -> %d)."
                        % (len(refs), len(deduped)))

    m = dict(DEFAULTS)
    for k in DEFAULTS:
        if not _blank(params.get(k)):
            m[k] = params[k]
    m["workspace_name"] = str(ws).strip()
    m["assembly_refs"] = deduped

    m["species"] = str(m["species"]).strip()
    m["custom_library_ref"] = str(m["custom_library_ref"]).strip()
    if not m["species"] and not m["custom_library_ref"]:
        raise ValueError("Give a species/taxon (e.g. 'Drosophila melanogaster', "
                         "'Fungi', 'Viridiplantae') to search Dfam, and/or a custom "
                         "repeat library Assembly. Neither was provided.")
    if m["species"] and not SPECIES_RE.match(m["species"]):
        raise ValueError("Species/taxon %r contains unsupported characters. Use an "
                         "NCBI Taxonomy scientific name such as 'Arabidopsis "
                         "thaliana' or a clade such as 'Fungi'." % m["species"])

    for key, allowed in (("mask_type", MASK_TYPES), ("mask_scope", MASK_SCOPES),
                         ("speed", SPEEDS)):
        m[key] = str(m[key]).strip()
        if m[key] not in allowed:
            raise ValueError("Parameter '%s' must be one of %s, got %r"
                             % (key, ", ".join(allowed), m[key]))

    for key in ("norna", "uncurated", "alignments", "save_masked_assembly"):
        m[key] = 1 if _as_bool(m[key]) else 0

    if m["div"] is not None:
        m["div"] = _as_float(m["div"], "div")
        if not 0 < m["div"] < 100:
            raise ValueError("Parameter 'div' (max %% divergence) must be between 0 "
                             "and 100, got %s" % m["div"])
    m["cutoff"] = _as_int(m["cutoff"], "cutoff")
    if not 100 <= m["cutoff"] <= 2000:
        raise ValueError("Parameter 'cutoff' must be between 100 and 2000, got %d"
                         % m["cutoff"])

    m["output_suffix"] = str(m["output_suffix"]).strip()
    if not SUFFIX_RE.match(m["output_suffix"]):
        raise ValueError("Parameter 'output_suffix' may only contain letters, digits "
                         "and . _ - | (max 40 characters), got %r" % m["output_suffix"])
    if m["save_masked_assembly"] and m["output_suffix"] == "":
        m["output_suffix"] = "_masked"
        warnings.append("Empty output suffix replaced by '_masked' so the input "
                        "objects are never overwritten.")

    # Library mode, and parameters that cannot take effect in it (lessons §9).
    if m["species"] and m["custom_library_ref"]:
        m["library_mode"] = "combined"
    elif m["custom_library_ref"]:
        m["library_mode"] = "custom"
    else:
        m["library_mode"] = "dfam"
    if m["library_mode"] == "dfam" and m["cutoff"] != DEFAULT_CUTOFF:
        warnings.append("'cutoff' only applies to custom-library searches; it was "
                        "not passed (Dfam search uses its own per-family thresholds).")
        m["cutoff"] = DEFAULT_CUTOFF
    if m["library_mode"] == "custom" and m["uncurated"]:
        warnings.append("'Include uncurated Dfam families' has no effect without a "
                        "species/taxon; it was not passed.")
        m["uncurated"] = 0
    if m["mask_scope"] == "simple_only":
        if m["custom_library_ref"]:
            warnings.append("Mask scope 'simple repeats / low complexity only' "
                            "skips interspersed repeats, so the custom library is "
                            "not searched.")
        if m["div"] is not None:
            warnings.append("'div' only filters interspersed repeats; it has no "
                            "effect with mask scope 'simple only'.")
    return m, warnings


def build_rm_command(rm_exe, input_fasta, out_dir, m, parallel, lib_path=None):
    """
    RepeatMasker [options] <seqfile>.  The sequence file must be LAST and every
    path absolute (RepeatMasker runs with cwd = job dir).
    """
    cmd = [rm_exe, "-engine", "rmblast", "-pa", str(parallel),
           "-dir", os.path.abspath(out_dir), "-gff", "-excln"]
    if lib_path:
        cmd += ["-lib", os.path.abspath(lib_path)]
        if m["cutoff"] != DEFAULT_CUTOFF:
            cmd += ["-cutoff", str(m["cutoff"])]
    else:
        cmd += ["-species", m["species"]]
        if m["uncurated"]:
            cmd += ["-uncurated"]
    if m["mask_type"] == "soft":
        cmd += ["-xsmall"]
    elif m["mask_type"] == "hard_X":
        cmd += ["-x"]
    if m["mask_scope"] == "interspersed_only":
        cmd += ["-nolow"]
    elif m["mask_scope"] == "simple_only":
        cmd += ["-noint"]
    cmd += {"quick": ["-q"], "rush": ["-qq"], "slow": ["-s"]}.get(m["speed"], [])
    if m["norna"]:
        cmd += ["-norna"]
    if m["div"] is not None:
        cmd += ["-div", ("%g" % m["div"])]
    if m["alignments"]:
        cmd += ["-a"]
    cmd.append(os.path.abspath(input_fasta))
    return cmd


def default_parallel():
    try:
        ncpu = len(os.sched_getaffinity(0))
    except AttributeError:
        ncpu = os.cpu_count() or 1
    return max(1, ncpu // 4)  # RMBlast uses 4 cores per -pa job


class RepeatMaskerRunner:

    def __init__(self, config, callback_url=None, ctx=None, clients=None):
        self.config = config
        self.callback_url = callback_url
        self.ctx = ctx
        self.scratch = os.path.abspath(config["scratch"])
        self._clients = clients or {}
        self.commands = []
        self.results = []
        self.log = logging.getLogger("kb_repeatmasker")

    # ------------------------------------------------------------------ clients
    def _client(self, name):
        if name not in self._clients:
            token = self.ctx["token"] if self.ctx else None
            if name == "au":
                from installed_clients.AssemblyUtilClient import AssemblyUtil
                self._clients[name] = AssemblyUtil(self.callback_url, token=token)
            elif name == "report":
                from installed_clients.KBaseReportClient import KBaseReport
                self._clients[name] = KBaseReport(self.callback_url, token=token)
        return self._clients[name]

    # --------------------------------------------------------------- subprocess
    def _env(self):
        env = dict(os.environ)
        env["PATH"] = ENV_BIN + os.pathsep + env.get("PATH", "")
        env["LC_ALL"] = "C"
        env["LANG"] = "C"
        # (RepeatMasker finds FamDB via famdb.conf, written in the Dockerfile.)
        home = os.path.join(self.scratch, "rm_home")
        os.makedirs(home, exist_ok=True)
        env["HOME"] = home  # fallback cache (~/.RepeatMaskerCache) must be writable
        return env

    def run_cmd(self, cmd, log_path, cwd, record=True, ok_markers=(), stdout_path=None):
        """Run with isolated-env PATH, log everything, raise with a log tail."""
        if record:
            self.commands.append(" ".join(shlex.quote(c) for c in cmd))
        tail = collections.deque(maxlen=60)
        captured = []
        with open(log_path, "a") as lf:
            lf.write("$ %s\n" % " ".join(shlex.quote(c) for c in cmd))
            lf.flush()
            if stdout_path:
                # stdout is data (e.g. a FASTA export); stderr goes to the log.
                with open(stdout_path, "w") as out_fh:
                    rc = subprocess.call(cmd, cwd=cwd, env=self._env(),
                                         stdout=out_fh, stderr=lf)
            else:
                p = subprocess.Popen(cmd, cwd=cwd, env=self._env(),
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     universal_newlines=True, bufsize=1)
                for line in p.stdout:
                    lf.write(line)
                    captured.append(line)
                    tail.append(line.rstrip("\n"))
                rc = p.wait()
        if stdout_path and rc != 0:
            with open(log_path) as lf:
                tail.extend(lf.read().splitlines()[-40:])
        text = "".join(captured)
        if rc != 0 and not any(mk in text for mk in ok_markers):
            raise RuntimeError("%s failed (exit %d):\n%s"
                               % (os.path.basename(cmd[0]), rc, "\n".join(tail)))
        return rc, text

    # -------------------------------------------------------------- tool facts
    def tool_versions(self, job_dir, log):
        notes, warnings, facts = [], [], {}
        _rc, text = self.run_cmd([os.path.join(ENV_BIN, "RepeatMasker"), "-v"],
                                 log, job_dir, record=False)
        facts["repeatmasker_version"] = P.parse_rm_version(text) or "unknown"
        if facts["repeatmasker_version"] != EXPECTED_RM_VERSION:
            warnings.append("RepeatMasker version is %s; this module was validated "
                            "with %s and Dfam %s (FamDB format 3 needs >= 4.2.4)."
                            % (facts["repeatmasker_version"], EXPECTED_RM_VERSION,
                               EXPECTED_DFAM_VERSION))
        return facts, notes, warnings

    def famdb_ready(self):
        return (os.path.isdir(FAMDB_DATA_DIR)
                and bool(glob.glob(os.path.join(FAMDB_DATA_DIR, "*.0.h5"))))

    def famdb(self, args, job_dir, log, record=True, stdout_path=None):
        cmd = [os.path.join(ENV_BIN, "famdb.py"), "-i", FAMDB_DATA_DIR] + args
        return self.run_cmd(cmd, log, job_dir, record=record, stdout_path=stdout_path)

    def dfam_facts(self, job_dir, log):
        _rc, text = self.famdb(["info"], job_dir, log, record=False)
        info = P.parse_famdb_info(text)
        files = sorted(os.path.basename(p)
                       for p in glob.glob(os.path.join(FAMDB_DATA_DIR, "*.h5")))
        info["files"] = files
        manifest = os.path.join(REFDATA_ROOT, "famdb_manifest.tsv")
        info["manifest"] = open(manifest).read() if os.path.exists(manifest) else ""
        return info

    def check_species(self, m, job_dir, log):
        """Resolve the taxon in FamDB before RepeatMasker sees it."""
        notes, warnings = [], []
        flags = [] if m["uncurated"] else ["-c"]
        _rc, text = self.famdb(["lineage", "-ad"] + flags + ["--format", "totals",
                               m["species"]], job_dir, log)
        totals = P.parse_lineage_totals(text)
        if totals is None:
            hint = "\n".join(text.strip().splitlines()[-15:])
            raise ValueError(
                "Species/taxon '%s' was not found unambiguously in Dfam (FamDB). "
                "Use an NCBI Taxonomy scientific name, e.g. 'Drosophila "
                "melanogaster', 'Arabidopsis thaliana', 'Saccharomycetes', 'Fungi'. "
                "FamDB said:\n%s" % (m["species"], hint or "(no output)"))
        present = (totals["ancestral_present"], totals["lineage_specific_present"])
        notes.append("Dfam families for '%s' (%s): %d ancestral + %d lineage-specific"
                     " in the release%s." % (
                         m["species"],
                         "curated + uncurated" if m["uncurated"] else "curated only",
                         totals["ancestral"], totals["lineage_specific"],
                         "" if present[0] is None else
                         "; %d + %d installed here" % present))
        n_avail = (sum(present) if present[0] is not None
                   else totals["ancestral"] + totals["lineage_specific"])
        if n_avail == 0:
            warnings.append("No Dfam families are available for '%s'. Only simple "
                            "repeats and low-complexity DNA will be masked. Consider "
                            "a custom repeat library (e.g. from RepeatModeler), or "
                            "'Include uncurated families'." % m["species"])
        comps = ["--component", "cc"] + (["--component", "uc"] if m["uncurated"] else [])
        _rc, ctext = self.famdb(["check"] + comps + [m["species"]], job_dir, log,
                                record=False)
        missing = P.parse_famdb_check(ctext)
        if missing:
            warnings.append(
                "Dfam partitions needed for '%s' are not installed in this module's "
                "reference data, so some families were not searched: %s. Add them "
                "(as .h5.gz names) to data/refdata_famdb.txt, bump data-version in "
                "kbase.yml and re-register." % (m["species"], ", ".join(missing)))
        return totals, notes, warnings

    # ----------------------------------------------------------------- library
    def prepare_library(self, m, job_dir, log):
        """Returns (lib_path or None, notes, warnings)."""
        notes, warnings = [], []
        if m["library_mode"] == "dfam":
            return None, notes, warnings
        lib_dir = os.path.join(job_dir, "library")
        os.makedirs(lib_dir, exist_ok=True)
        got = self._client("au").get_assembly_as_fasta({"ref": m["custom_library_ref"]})
        custom = os.path.join(lib_dir, "custom_library.fa")
        n_custom = n_class = 0
        with open(custom, "w") as out:
            for sid, _d, seq in P.iter_fasta(got["path"]):
                n_custom += 1
                n_class += "#" in sid
                P.write_fasta_record(out, sid, seq.upper())
        if n_custom == 0:
            raise ValueError("The custom repeat library %s contains no sequences"
                             % m["custom_library_ref"])
        lib_name = got.get("assembly_name") or m["custom_library_ref"]
        notes.append("Custom repeat library '%s': %d sequences." % (lib_name, n_custom))
        if n_class < n_custom:
            warnings.append("%d of %d custom library sequences have no '#class/family' "
                            "in their ID (RepeatMasker format 'name#LTR/Gypsy'); they "
                            "are reported as Unknown." % (n_custom - n_class, n_custom))
        if m["library_mode"] == "custom":
            return custom, notes, warnings
        # combined: Dfam families for the taxon + custom library (repeatmasker.help §1.1)
        dfam_lib = os.path.join(lib_dir, "dfam_%s.fa" % safe_name(m["species"]))
        args = ["families", "-f", "fasta_name", "--include-class-in-name", "-ad"]
        if not m["uncurated"]:
            args.append("--curated")
        self.famdb(args + [m["species"]], job_dir, log, stdout_path=dfam_lib)
        n_dfam = P.count_fasta_records(dfam_lib)
        combined = os.path.join(lib_dir, "combined_library.fa")
        with open(combined, "w") as out:
            for path in (dfam_lib, custom):
                with open(path) as fh:
                    shutil.copyfileobj(fh, out)
        notes.append("Combined library: %d Dfam families for '%s' + %d custom "
                     "sequences, searched with -lib (RepeatMasker's staged mammalian "
                     "search is not used in this mode)." % (n_dfam, m["species"],
                                                            n_custom))
        if n_dfam == 0:
            warnings.append("FamDB exported 0 families for '%s'; only the custom "
                            "library was searched." % m["species"])
        return combined, notes, warnings

    # -------------------------------------------------------------------- run
    def run(self, params):
        m, warnings = validate_params(params)
        notes = []
        job_dir = os.path.join(self.scratch, "rm_" + uuid.uuid4().hex[:12])
        os.makedirs(job_dir)
        log = os.path.join(job_dir, "repeatmasker_setup.log")
        self.commands = []
        self.results = []

        # RepeatMasker 4.2.4 runs `famdb.py info` on EVERY run, even with -lib only,
        # so FamDB data is required in all library modes.
        if not self.famdb_ready():
            raise RuntimeError(
                "Dfam reference data was not found in %s. The module's reference "
                "data (entrypoint.sh init) has not been installed; ask the module "
                "owner to check the catalog registration log." % FAMDB_DATA_DIR)

        facts, n2, w2 = self.tool_versions(job_dir, log)
        notes += n2
        warnings += w2
        dfam = {}
        if m["library_mode"] in ("dfam", "combined"):
            dfam = self.dfam_facts(job_dir, log)
            if dfam.get("db_version") and dfam["db_version"] != EXPECTED_DFAM_VERSION:
                warnings.append("Installed Dfam is %s; this module expects %s."
                                % (dfam["db_version"], EXPECTED_DFAM_VERSION))
            _t, n3, w3 = self.check_species(m, job_dir, log)
            notes += n3
            warnings += w3
            notes.append("Dfam source: module reference data (%s): %s."
                         % (FAMDB_DATA_DIR, ", ".join(dfam.get("files", []))))
        lib_path, n4, w4 = self.prepare_library(m, job_dir, log)
        notes += n4
        warnings += w4

        parallel = default_parallel()
        notes.append("RepeatMasker ran with -pa %d (RMBlast uses 4 threads per job)."
                     % parallel)
        notes.append("Percentages exclude runs of >= 20 N/X (-excln), so scaffold gaps"
                     " do not dilute repeat content.")
        if m["mask_scope"] == "interspersed_only":
            notes.append("-nolow: simple repeats/low complexity were not searched at "
                         "all, so short simple stretches inside TEs may be annotated "
                         "as TE fragments (RepeatMasker documentation).")

        out_refs = []
        used_names = set()
        for idx, ref in enumerate(m["assembly_refs"], start=1):
            res = self.run_one(idx, ref, m, lib_path, parallel, job_dir, used_names)
            self.results.append(res)
            if res.get("masked_assembly_ref"):
                out_refs.append(res["masked_assembly_ref"])

        report = self.make_report(m, facts, dfam, notes, warnings, job_dir)
        compact = [{k: r.get(k) for k in (
            "input_ref", "input_name", "masked_assembly_ref", "masked_assembly_name",
            "contigs", "total_bp", "masked_bp", "masked_pct", "interspersed_bp",
            "interspersed_pct", "simple_lowcomp_bp", "n_annotations")}
            for r in self.results]
        for c, r in zip(compact, self.results):
            c["warnings"] = len(r["warnings"])
            c["categories_bp"] = {k: v["bp"] for k, v in r["categories"].items()}
        return {"report_name": report["name"], "report_ref": report["ref"],
                "masked_assembly_refs": out_refs, "repeatmasker_results": compact}

    def run_one(self, idx, ref, m, lib_path, parallel, job_dir, used_names):
        warnings = []
        got = self._client("au").get_assembly_as_fasta({"ref": ref})
        name = got.get("assembly_name") or ref.replace("/", "_")
        tag = "asm%02d_%s" % (idx, safe_name(name, 40))
        work = os.path.join(job_dir, tag)
        raw = os.path.join(work, "raw")
        os.makedirs(raw)
        log = os.path.join(work, "repeatmasker.log")

        query = os.path.join(raw, "query.fa")  # safe name: no spaces / specials
        mapping = P.rename_fasta(os.path.abspath(got["path"]), query)
        id_map = {x["short_id"]: x["original_id"] for x in mapping}
        P.write_id_map(mapping, os.path.join(work, "contig_id_map.tsv"))
        total_bp, gc = P.gc_and_length(query)

        cmd = build_rm_command(os.path.join(ENV_BIN, "RepeatMasker"), query, raw, m,
                               parallel, lib_path)
        _rc, text = self.run_cmd(cmd, log, raw, ok_markers=P.NO_REPEATS_MARKERS)
        no_repeats_msg = any(mk in text for mk in P.NO_REPEATS_MARKERS)

        base = os.path.join(raw, "query.fa")
        out_path, tbl_path = base + ".out", base + ".tbl"
        masked_raw = base + ".masked"
        # RepeatMasker 4.2.4 writes <seq>.out.gff (GFF3); older releases <seq>.gff
        gff_raw = next((g for g in (base + ".out.gff", base + ".gff")
                        if os.path.exists(g)), base + ".out.gff")
        if not os.path.exists(out_path):
            raise RuntimeError("RepeatMasker produced no .out file for %s; see the "
                               "log:\n%s" % (name, "\n".join(text.splitlines()[-30:])))
        rows, no_rep_out = P.parse_out(out_path, id_map)
        safe = safe_name(name)
        out_tr = os.path.join(work, safe + ".out")
        P.translate_out_file(out_path, out_tr, id_map)
        if os.path.exists(gff_raw):
            P.translate_gff_file(gff_raw, os.path.join(work, safe + ".gff"), id_map)
        if os.path.exists(tbl_path):
            shutil.copy(tbl_path, os.path.join(work, safe + ".tbl"))

        masked_fa = os.path.join(work, safe + ".masked.fa")
        if os.path.exists(masked_raw):
            per_contig, problems = P.translate_masked_fasta(
                got["path"], masked_raw, masked_fa, mapping, m["mask_type"])
            warnings += ["Integrity: " + p for p in problems]
        else:
            if rows:
                raise RuntimeError("RepeatMasker reported %d annotations for %s but "
                                   "wrote no .masked file" % (len(rows), name))
            per_contig = P.copy_unmasked(got["path"], masked_fa)
        with open(os.path.join(work, "per_contig_masking.tsv"), "w") as fh:
            fh.write("contig_id\tlength\tmasked_bp\tmasked_pct\n")
            for c in per_contig:
                fh.write("%s\t%d\t%d\t%.2f\n" % (c["contig_id"], c["length"],
                         c["masked_bp"], 100.0 * c["masked_bp"] / max(1, c["length"])))

        masked_bp = sum(c["masked_bp"] for c in per_contig)
        summ = P.summarize_rows(rows)
        tbl = P.parse_tbl(tbl_path) if os.path.exists(tbl_path) else {}
        if "bases_masked" in tbl and masked_bp:
            diff = abs(tbl["bases_masked"] - masked_bp) / float(masked_bp)
            if diff > 0.02:
                warnings.append("Masked bases in the masked FASTA (%d) and the .tbl "
                                "summary (%d) differ by %.1f%%."
                                % (masked_bp, tbl["bases_masked"], 100 * diff))
        if not rows:
            warnings.append("No repeats were annotated in %s%s." % (
                name, " (RepeatMasker: no repetitive sequences detected)"
                if (no_repeats_msg or no_rep_out) else ""))
        n_contigs = len(mapping)
        if len(per_contig) != n_contigs:
            warnings.append("Input had %d contigs but the masked output has %d."
                            % (n_contigs, len(per_contig)))
        pct = 100.0 * masked_bp / max(1, total_bp)
        if m["library_mode"] != "dfam" and summ["interspersed_bp"] > 0.8 * total_bp:
            warnings.append("Over 80%% of %s is annotated as interspersed repeats; "
                            "check the custom library for genes or low-complexity "
                            "sequence." % name)
        alert = base + ".alert"
        if os.path.exists(alert) and os.path.getsize(alert):
            warnings.append("RepeatMasker flagged possible E. coli insertion "
                            "elements/contamination (see the .alert file in the zip).")
            shutil.copy(alert, os.path.join(work, safe + ".alert"))

        masked_ref = masked_obj = None
        if m["save_masked_assembly"]:
            masked_obj = safe_name(name, 200) + m["output_suffix"]
            k = 2
            while masked_obj in used_names:
                masked_obj = "%s%s_%d" % (safe_name(name, 200), m["output_suffix"], k)
                k += 1
            used_names.add(masked_obj)
            to_save = os.path.join(self.scratch, "%s_%s.fa" % (tag, uuid.uuid4().hex[:6]))
            shutil.copy(masked_fa, to_save)
            masked_ref = self._client("au").save_assembly_from_fasta({
                "file": {"path": to_save}, "workspace_name": m["workspace_name"],
                "assembly_name": masked_obj})

        # Raw files with temporary rmq IDs are superseded by the translated copies.
        for raw_f in (masked_raw, out_path, tbl_path, gff_raw):
            if os.path.exists(raw_f):
                os.remove(raw_f)
        os.remove(query)  # input copy; contig_id_map.tsv documents the renaming
        zip_path = shutil.make_archive(os.path.join(job_dir, tag), "zip",
                                       root_dir=job_dir, base_dir=tag)
        return {
            "index": idx, "tag": tag, "input_ref": ref, "input_name": name,
            "masked_assembly_ref": masked_ref, "masked_assembly_name": masked_obj,
            "contigs": n_contigs, "total_bp": total_bp, "gc_pct": round(gc, 2),
            "masked_bp": masked_bp, "masked_pct": round(pct, 2),
            "interspersed_bp": summ["interspersed_bp"],
            "interspersed_pct": round(100.0 * summ["interspersed_bp"] / max(1, total_bp), 2),
            "simple_lowcomp_bp": summ["simple_lowcomp_bp"],
            "n_annotations": len(rows), "categories": summ["categories"],
            "top_families": summ["top_families"], "n_families": summ["n_families"],
            "tbl": tbl, "rows": rows, "per_contig": per_contig,
            "warnings": warnings, "zip": zip_path, "masked_fasta": masked_fa,
            "work_dir": work,
        }

    # ------------------------------------------------------------------ report
    def write_summary_tsv(self, path):
        cats = list(P.CATEGORY_ORDER)
        head = ["input_name", "input_ref", "masked_assembly_name", "masked_assembly_ref",
                "contigs", "total_bp", "gc_pct", "masked_bp", "masked_pct",
                "interspersed_bp", "interspersed_pct", "simple_lowcomp_bp",
                "annotations"] + [c.replace(" ", "_") + "_bp" for c in cats]
        with open(path, "w") as fh:
            fh.write("\t".join(head) + "\n")
            for r in self.results:
                vals = [r["input_name"], r["input_ref"], r["masked_assembly_name"] or "",
                        r["masked_assembly_ref"] or "", r["contigs"], r["total_bp"],
                        r["gc_pct"], r["masked_bp"], r["masked_pct"],
                        r["interspersed_bp"], r["interspersed_pct"],
                        r["simple_lowcomp_bp"], r["n_annotations"]]
                vals += [r["categories"].get(c, {"bp": 0})["bp"] for c in cats]
                fh.write("\t".join(str(v) for v in vals) + "\n")

    def make_report(self, m, facts, dfam, notes, warnings, job_dir):
        report_dir = os.path.join(job_dir, "report")
        os.makedirs(report_dir)
        tsv = os.path.join(job_dir, "repeatmasker_summary.tsv")
        self.write_summary_tsv(tsv)
        with open(os.path.join(report_dir, "index.html"), "w") as fh:
            fh.write(render_html(m, facts, dfam, notes, warnings, self.results,
                                 self.commands))
        file_links = [{"path": tsv, "name": "repeatmasker_summary.tsv",
                       "description": "Per-assembly repeat content summary (TSV)"}]
        for r in self.results:
            file_links.append({"path": r["zip"], "name": os.path.basename(r["zip"]),
                               "description": "Full RepeatMasker output for %s "
                               "(.out, .gff, .tbl, masked FASTA, logs)" % r["input_name"]})
        objs = [{"ref": r["masked_assembly_ref"],
                 "description": "%s-masked assembly of %s" % (
                     m["mask_type"].replace("_", "-"), r["input_name"])}
                for r in self.results if r["masked_assembly_ref"]]
        total_masked = sum(r["masked_bp"] for r in self.results)
        msg = "RepeatMasker %s: %d assembl%s, %s bp masked." % (
            facts.get("repeatmasker_version", ""), len(self.results),
            "y" if len(self.results) == 1 else "ies", format(total_masked, ","))
        params = {
            "message": msg, "objects_created": objs,
            "direct_html_link_index": 0,
            "html_links": [{"path": report_dir, "name": "index.html",
                            "description": "RepeatMasker report"}],
            "file_links": file_links,
            "report_object_name": "kb_repeatmasker_report_" + uuid.uuid4().hex[:8],
            "workspace_name": m["workspace_name"],
            "html_window_height": 800,
        }
        return self._client("report").create_extended_report(params)


# ---------------------------------------------------------------------- HTML
_COLORS = {"SINE": "#4e79a7", "LINE": "#f28e2b", "LTR": "#e15759", "DNA": "#76b7b2",
           "RC": "#59a14f", "Retroposon": "#edc948", "Unknown": "#b07aa1",
           "Other": "#ff9da7", "Satellite": "#9c755f", "Small RNA": "#bab0ac",
           "Simple_repeat": "#8cd17d", "Low_complexity": "#d4a6c8"}


def _e(x):
    return html.escape(str(x))


def _fmt(n):
    return format(int(n), ",")


def render_html(m, facts, dfam, notes, warnings, results, commands):
    all_w = list(warnings)
    for r in results:
        all_w += ["%s: %s" % (r["input_name"], w) for w in r["warnings"]]
    css = """
    body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:18px;
         color:#222;max-width:1150px}
    h1{font-size:20px;margin:0 0 4px} h2{font-size:16px;margin:22px 0 8px;
       border-bottom:1px solid #ddd;padding-bottom:4px} h3{font-size:14px;margin:14px 0 6px}
    table{border-collapse:collapse;font-size:12.5px;margin-bottom:8px}
    th,td{border:1px solid #ddd;padding:4px 8px;text-align:right}
    th{background:#f4f6f8} td.l,th.l{text-align:left}
    .warn{background:#fff4e5;border-left:4px solid #f0a020;padding:8px 12px}
    .note{background:#eef5fb;border-left:4px solid #4e79a7;padding:8px 12px}
    .bar{display:flex;height:18px;width:100%;max-width:700px;border:1px solid #ccc;
         background:#fafafa} .bar div{height:100%}
    .legend span{display:inline-block;margin-right:12px;font-size:12px}
    .legend i{display:inline-block;width:10px;height:10px;margin-right:4px}
    pre{background:#f6f8fa;padding:8px;font-size:11.5px;white-space:pre-wrap;
        word-break:break-all} .muted{color:#666;font-size:12px}
    """
    h = ["<!DOCTYPE html><html><head><meta charset='utf-8'><title>RepeatMasker "
         "report</title><style>%s</style></head><body>" % css]
    lib_desc = {"dfam": "Dfam (FamDB) families for species/taxon '%s'" % m["species"],
                "custom": "custom repeat library only",
                "combined": "Dfam families for '%s' + custom library" % m["species"]}
    sub = [lib_desc[m["library_mode"]]]
    if dfam.get("db_version"):
        sub.append("Dfam %s%s" % (dfam["db_version"], " (%s)" % dfam["db_date"]
                                  if dfam.get("db_date") else ""))
    sub += ["masking: " + m["mask_type"], "scope: " + m["mask_scope"],
            "speed: " + m["speed"],
            "run " + datetime.datetime.now(datetime.timezone.utc)
            .strftime("%Y-%m-%d %H:%M UTC")]
    h.append("<h1>RepeatMasker %s</h1><div class='muted'>%s</div>"
             % (_e(facts.get("repeatmasker_version", "")),
                " &middot; ".join(_e(x) for x in sub)))
    if all_w:
        h.append("<h2>Warnings</h2><div class='warn'><ul>%s</ul></div>"
                 % "".join("<li>%s</li>" % _e(w) for w in all_w))

    h.append("<h2>Summary</h2><table><tr><th class='l'>Assembly</th><th>Contigs</th>"
             "<th>Length (bp)</th><th>GC %</th><th>Masked (bp)</th><th>Masked %</th>"
             "<th>Interspersed %</th><th>Simple + low-complexity (bp)</th>"
             "<th>Annotations</th><th class='l'>Masked assembly</th></tr>")
    for r in results:
        h.append("<tr><td class='l'>%s</td><td>%s</td><td>%s</td><td>%.2f</td><td>%s</td>"
                 "<td>%.2f</td><td>%.2f</td><td>%s</td><td>%s</td><td class='l'>%s</td></tr>"
                 % (_e(r["input_name"]), _fmt(r["contigs"]), _fmt(r["total_bp"]),
                    r["gc_pct"], _fmt(r["masked_bp"]), r["masked_pct"],
                    r["interspersed_pct"], _fmt(r["simple_lowcomp_bp"]),
                    _fmt(r["n_annotations"]), _e(r["masked_assembly_name"] or "not saved")))
    h.append("</table>")

    legend = "".join("<span><i style='background:%s'></i>%s</span>" % (c, _e(k))
                     for k, c in _COLORS.items())
    for r in results:
        h.append("<h2>%s</h2>" % _e(r["input_name"]))
        segs = []
        for k in P.CATEGORY_ORDER:
            bp = r["categories"].get(k, {"bp": 0})["bp"]
            if bp:
                w = 100.0 * bp / max(1, r["total_bp"])
                segs.append("<div title='%s: %.2f%%' style='width:%.3f%%;background:%s'>"
                            "</div>" % (_e(k), w, min(w, 100), _COLORS[k]))
        h.append("<div class='muted'>Share of the assembly annotated per class "
                 "(whole bar = 100%% of %s bp)</div><div class='bar'>%s</div>"
                 "<div class='legend'>%s</div>" % (_fmt(r["total_bp"]), "".join(segs),
                                                    legend))
        h.append("<h3>By class</h3><table><tr><th class='l'>Class</th><th>Elements</th>"
                 "<th>Annotated bp</th><th>% of assembly</th></tr>")
        for k in P.CATEGORY_ORDER:
            v = r["categories"].get(k)
            if v:
                h.append("<tr><td class='l'>%s</td><td>%s</td><td>%s</td><td>%.2f</td></tr>"
                         % (_e(k), _fmt(v["count"]), _fmt(v["bp"]),
                            100.0 * v["bp"] / max(1, r["total_bp"])))
        h.append("</table>")
        if r["top_families"]:
            h.append("<h3>Top families by annotated length (%d families in total)</h3>"
                     "<table><tr><th class='l'>Repeat</th><th class='l'>Class/family</th>"
                     "<th>Elements</th><th>bp</th></tr>" % r["n_families"])
            for f in r["top_families"]:
                h.append("<tr><td class='l'>%s</td><td class='l'>%s</td><td>%s</td>"
                         "<td>%s</td></tr>" % (_e(f["repeat"]), _e(f["class"]),
                                                _fmt(f["count"]), _fmt(f["bp"])))
            h.append("</table>")
        if r["tbl"]:
            h.append("<div class='muted'>RepeatMasker .tbl cross-check: bases masked "
                     "%s bp (%s %%).</div>" % (_fmt(r["tbl"].get("bases_masked", 0)),
                                                _e(r["tbl"].get("pct_masked", "?"))))

    if notes:
        h.append("<h2>Notes</h2><div class='note'><ul>%s</ul></div>"
                 % "".join("<li>%s</li>" % _e(n) for n in notes))
    shown = {k: v for k, v in m.items() if k not in ("workspace_name",)}
    h.append("<h2>Parameters</h2><pre>%s</pre>" % _e(json.dumps(shown, indent=1)))
    h.append("<h2>Commands executed</h2><pre>%s</pre>" % _e("\n\n".join(commands)))
    if dfam.get("manifest"):
        h.append("<h2>Reference data</h2><pre>%s</pre>" % _e(dfam["manifest"]))
    h.append("<h2>Citation</h2><ul>%s</ul>" % "".join("<li>%s</li>" % _e(c)
                                                        for c in CITATIONS))
    h.append("<p class='muted'>Masking: soft = repeats in lowercase (recommended for "
             "gene prediction, e.g. BRAKER/AUGUSTUS); hard = repeats replaced by N or X."
             " Contig IDs were temporarily renamed (RepeatMasker limits IDs to 50 "
             "characters); all outputs use the original IDs (see contig_id_map.tsv)."
             "</p></body></html>")
    return "\n".join(h)
