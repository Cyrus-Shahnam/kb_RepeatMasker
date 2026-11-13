#!/bin/bash
# Read-only diagnostic block (lessons §0). Run from the module root and paste the
# whole output back when asking for help. Changes nothing.
cd "$(dirname "$0")/.." || exit 1
echo "== pwd (must contain kbase.yml)"; pwd; ls kbase.yml
echo "== UI dirs (want only run_repeatmasker)"; ls ui/narrative/methods/
echo "== funcdef"; grep -n funcdef kb_repeatmasker.spec
echo "== Impl delegates"; grep -c RepeatMaskerRunner lib/kb_repeatmasker/kb_repeatmaskerImpl.py
echo "== Java flags (want none)"; grep -n -- '--java' Makefile || echo none
echo "== clients (need AssemblyUtil + KBaseReport)"; ls lib/installed_clients/
echo "== empty files (want none)"
find . -type f -empty -not -path './test_local/*' -not -path './.git/*' -not -name '__init__.py'
echo "== first lines"
for f in kbase.yml kb_repeatmasker.spec Dockerfile env-repeatmasker.yml scripts/entrypoint.sh \
         scripts/download_dfam_data.sh lib/kb_repeatmasker/*.py test/*.py tools/*.py \
         ui/narrative/methods/run_repeatmasker/spec.json; do
  printf '%-52s %5s  %s\n' "$f" "$(wc -l < "$f" 2>/dev/null)" "$(head -1 "$f" 2>/dev/null | cut -c1-40)"
done
echo "== heredoc wrappers pasted into files (want none)"
grep -lE '^(cat > |EOF$)' kbase.yml kb_repeatmasker.spec Dockerfile scripts/*.sh lib/kb_repeatmasker/*.py test/*.py 2>/dev/null || echo none
echo "== exec bits"; ls -l scripts/*.sh tools/fake_bin/* | awk '{print $1, $NF}'
echo "== test.cfg"; ls test_local/test.cfg 2>/dev/null && \
  awk -F= '/^test_token=/{print "token length:", length($2)}' test_local/test.cfg && \
  echo "escape chars: $(grep -c $'\e' test_local/test.cfg)" && grep -E '^kbase_endpoint' test_local/test.cfg
echo "== token not tracked by git"; git check-ignore test_local/test.cfg 2>/dev/null || echo "WARNING: not ignored (or not a git repo)"
echo "== name chain"; python3 tools/check_names.py
