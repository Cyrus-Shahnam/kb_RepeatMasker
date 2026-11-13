#!/usr/bin/env python3
"""
Verify the hand-pasted kb_repeatmasker kit (lessons section 0, 15).
Run from the module root:  python3 tools/verify_kit.py
Checks every file: present, non-empty, expected line count (+0..+3 tolerated:
terminals add blank lines), expected first line (catches swapped files and
pasted heredoc wrappers), a parse check by type, executable bits, init
leftovers, and Java flags in the Makefile. Changes nothing.
"""
import ast
import json
import os
import subprocess
import sys

EXPECTED = {
 ".dockerignore": [
  5,
  [
   "test_local",
   "sdk.cfg",
   ".git"
  ]
 ],
 "Dockerfile": [
  70,
  [
   "FROM kbase/sdkpython:3.8.10",
   "LABEL maintainer=\"REPLACE_WITH_YOUR_NAME_OR_EMAIL\"",
   "USER root"
  ]
 ],
 "README.md": [
  82,
  [
   "# kb_repeatmasker",
   "",
   "KBase SDK module wrapping **RepeatMasker 4.2.4** (RMBlast en"
  ]
 ],
 "data/refdata_famdb.txt": [
  20,
  [
   "# Dfam FamDB files installed as KBase reference data (/data/",
   "# Read by scripts/download_dfam_data.sh during `entrypoint.s",
   "#"
  ]
 ],
 "env-repeatmasker.yml": [
  17,
  [
   "name: repeatmasker",
   "channels:",
   "  - conda-forge"
  ]
 ],
 "kb_repeatmasker.spec": [
  69,
  [
   "/*",
   "A KBase module: kb_repeatmasker",
   ""
  ]
 ],
 "kbase.yml": [
  24,
  [
   "#",
   "# define display information",
   "#"
  ]
 ],
 "lib/kb_repeatmasker/kb_repeatmaskerImpl.py": [
  74,
  [
   "# -*- coding: utf-8 -*-",
   "#BEGIN_HEADER",
   "import logging"
  ]
 ],
 "lib/kb_repeatmasker/repeatmasker_utils.py": [
  776,
  [
   "\"\"\"",
   "RepeatMaskerRunner: all tool logic for the kb_repeatmasker a",
   ""
  ]
 ],
 "lib/kb_repeatmasker/rm_parsers.py": [
  338,
  [
   "\"\"\"",
   "Pure-Python helpers for kb_repeatmasker: FASTA streaming, co",
   "RepeatMasker .out / .tbl / .gff parsing and masking statisti"
  ]
 ],
 "scripts/download_dfam_data.sh": [
  81,
  [
   "#!/bin/bash",
   "# Download the Dfam FamDB files listed in data/refdata_famdb",
   "# verify each against Dfam's .md5 sidecar, decompress, and p"
  ]
 ],
 "scripts/entrypoint.sh": [
  37,
  [
   "#!/bin/bash",
   "",
   "# sdkbase2 ships user-env.sh; sdkpython does not."
  ]
 ],
 "test/kb_repeatmasker_server_test.py": [
  339,
  [
   "# -*- coding: utf-8 -*-",
   "\"\"\"",
   "kb_repeatmasker tests."
  ]
 ],
 "test/make_fixtures.py": [
  62,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "Seeded known-truth fixtures for kb_repeatmasker (lessons \u00a78)"
  ]
 ],
 "tools/check_names.py": [
  59,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "Name-chain check (lessons \u00a72, \u00a717). Run from the module root"
  ]
 ],
 "tools/diagnose.sh": [
  26,
  [
   "#!/bin/bash",
   "# Read-only diagnostic block (lessons \u00a70). Run from the modu",
   "# whole output back when asking for help. Changes nothing."
  ]
 ],
 "tools/dry_run.py": [
  222,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "Dry-run the whole Runner OUTSIDE Docker with a fake RepeatMa"
  ]
 ],
 "tools/fake_bin/RepeatMasker": [
  168,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "FAKE RepeatMasker for dry-running the kb_repeatmasker Runner"
  ]
 ],
 "tools/fake_bin/famdb.py": [
  49,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "FAKE famdb.py for dry runs (see tools/fake_bin/RepeatMasker)"
  ]
 ],
 "tools/plant_repeats.py": [
  178,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "Plant known repeat copies into a genome, for known-truth tes"
  ]
 ],
 "tools/score_masking.py": [
  116,
  [
   "#!/usr/bin/env python3",
   "\"\"\"",
   "Score a masked FASTA against a plant_repeats.py manifest."
  ]
 ],
 "ui/narrative/methods/run_repeatmasker/display.yaml": [
  132,
  [
   "#",
   "# define display information",
   "#"
  ]
 ],
 "ui/narrative/methods/run_repeatmasker/spec.json": [
  328,
  [
   "{",
   "    \"ver\": \"0.1.0\",",
   "    \"authors\": ["
  ]
 ]
}
EXEC = ["scripts/entrypoint.sh", "scripts/download_dfam_data.sh", "tools/diagnose.sh",
        "tools/fake_bin/RepeatMasker", "tools/fake_bin/famdb.py"]
bad, warn = [], []


def line_ok(got, want):
    got = got.rstrip()[:60]
    if "REPLACE_WITH" in want:  # placeholder lines may already be filled in
        return got.startswith(want[:want.index("REPLACE_WITH")])
    return got == want


for path, (n, first) in sorted(EXPECTED.items()):
    if not os.path.exists(path):
        bad.append("MISSING  %s" % path)
        continue
    text = open(path, encoding="utf-8", errors="replace").read()
    if not text.strip():
        bad.append("EMPTY    %s" % path)
        continue
    got = text.count("\n")
    line1 = text.split("\n", 1)[0]
    if "\r" in text:
        bad.append("CRLF     %s (run: sed -i 's/\\r$//' %s)" % (path, path))
    if line1.startswith("cat >") or "\nEOF\n" in text or text.endswith("\nEOF\n"):
        bad.append("HEREDOC  %s (run: sed -i '1d;/^EOF$/,$d' %s)" % (path, path))
    elif not all(line_ok(g, e) for g, e in zip(text.split("\n")[:3] + [""] * 3, first)):
        bad.append("HEAD3    %s: first lines %r, want %r (wrong/swapped file?)"
                   % (path, text.split("\n")[:3], first))
    if path.endswith("Impl.py"):
        # `make` adds the spec's parameter docs to the Impl docstrings, so after a
        # compile the file is longer than the kit's. Check content instead.
        for mark in ("RepeatMaskerRunner(", "#BEGIN run_repeatmasker", "#END run_repeatmasker"):
            if mark not in text:
                bad.append("CONTENT  %s: %r missing" % (path, mark))
        if text.count("\nclass kb_repeatmasker") != 1:
            bad.append("CONTENT  %s: class defined %d times (duplicated paste?)"
                       % (path, text.count("\nclass kb_repeatmasker")))
        if got < n:
            bad.append("LINES    %s: %d lines, want at least %d" % (path, got, n))
    elif got < n or got > n + 3:
        bad.append("LINES    %s: %d lines, want %d (truncated or duplicated paste?)"
                   % (path, got, n))
    elif got != n and not path.endswith("Impl.py"):
        warn.append("lines    %s: %d vs %d (extra blank lines, harmless)" % (path, got, n))
    try:
        if path.endswith(".py") or path.endswith("fake_bin/RepeatMasker"):
            ast.parse(text)
        elif path.endswith(".json"):
            json.loads(text)
        elif path.endswith((".yml", ".yaml")):
            try:
                import yaml
                yaml.safe_load(text)
            except ImportError:
                warn.append("no PyYAML here; skipped YAML parse of %s" % path)
        elif path.endswith(".sh"):
            subprocess.check_call(["bash", "-n", path])
    except Exception as e:  # noqa
        bad.append("PARSE    %s: %s" % (path, str(e).splitlines()[0][:100]))
for path in EXEC:
    if os.path.exists(path) and not os.access(path, os.X_OK):
        bad.append("NOEXEC   %s (run: chmod +x %s)" % (path, path))
for d in os.listdir("ui/narrative/methods") if os.path.isdir("ui/narrative/methods") else []:
    if d != "run_repeatmasker":
        bad.append("LEFTOVER ui/narrative/methods/%s (delete init's example app)" % d)
for path in EXPECTED:  # IDE "fix imports" rewrites (PyCharm/VS Code)
    if path.endswith(".py") and os.path.exists(path):
        for n, line in enumerate(open(path, encoding="utf-8", errors="replace"), 1):
            if line.startswith(("from tools.", "from test.", "import tools", "from lib.")):
                bad.append("IMPORT   %s:%d %r (IDE rewrote an import; restore the kit line)"
                           % (path, n, line.strip()[:60]))
if os.path.exists("Makefile") and "--java" in open("Makefile").read():
    bad.append("JAVA     Makefile still has --java flags")
for need in ("AssemblyUtilClient.py", "KBaseReportClient.py"):
    if not os.path.exists(os.path.join("lib/installed_clients", need)):
        bad.append("CLIENT   lib/installed_clients/%s missing (kb-sdk install %s)"
                   % (need, need.replace("Client.py", "")))
if not os.path.exists("lib/kb_repeatmasker/authclient.py"):
    bad.append("AUTH     lib/kb_repeatmasker/authclient.py missing; the generated server "
               "imports it (run: cp lib/installed_clients/authclient.py lib/kb_repeatmasker/)")
for w in warn:
    print("note  " + w)
for b in bad:
    print("FAIL  " + b)
print("%d files checked, %d problems" % (len(EXPECTED), len(bad)))
sys.exit(1 if bad else 0)
