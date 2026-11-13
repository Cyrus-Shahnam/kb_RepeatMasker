#!/bin/bash
# Download the Dfam FamDB files listed in data/refdata_famdb.txt into <dest>/famdb,
# verify each against Dfam's .md5 sidecar, decompress, and prove FamDB can read
# them. Run by `entrypoint.sh init` (catalog registration, and kb-sdk test).
#   usage: download_dfam_data.sh <dest_dir>
set -euo pipefail

DEST="${1:?usage: download_dfam_data.sh <dest_dir>}"
MODULE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LIST="${REFDATA_LIST:-$MODULE_DIR/data/refdata_famdb.txt}"
ENV_BIN="${REPEATMASKER_ENV_BIN:-/opt/conda/envs/repeatmasker/bin}"
export PATH="$ENV_BIN:$PATH"            # this script only: famdb.py needs env python3
FAMDB="$DEST/famdb"
MD5DIR="$DEST/famdb_md5"                # kept OUT of the FamDB dir (only *.h5 there)
mkdir -p "$FAMDB" "$MD5DIR"

BASE_URL="${DFAM_BASE_URL:-$(awk '$1=="BASE_URL"{print $2}' "$LIST")}"
[ -n "$BASE_URL" ] || { echo "ERROR: no BASE_URL line in $LIST"; exit 1; }
FILES=$(awk '!/^#/ && NF && $1!="BASE_URL"{print $1}' "$LIST")
[ -n "$FILES" ] || { echo "ERROR: no files listed in $LIST"; exit 1; }

WORK="$(mktemp -d)"                     # tools may drop files in the cwd
cd "$WORK"
echo "Dfam source: $BASE_URL"
echo "Destination: $FAMDB"

for gz in $FILES; do
  stem="${gz%.gz}"
  if [ -s "$FAMDB/$stem" ] && [ -s "$MD5DIR/$gz.md5ok" ]; then
    echo "already present: $stem"; continue
  fi
  echo "== $gz"
  curl -fsSL --retry 10 --retry-delay 15 -o "$gz.md5" "$BASE_URL$gz.md5"
  expected=$(awk '{print tolower($1); exit}' "$gz.md5")
  [ ${#expected} -eq 32 ] || { echo "ERROR: bad md5 sidecar for $gz"; cat "$gz.md5"; exit 1; }
  actual=""
  for attempt in 1 2 3; do
    curl -fL --retry 10 --retry-delay 15 -C - -o "$gz" "$BASE_URL$gz" || true
    actual=$( [ -f "$gz" ] && md5sum "$gz" | awk '{print $1}' || true )
    [ "$actual" = "$expected" ] && break
    echo "md5 mismatch on attempt $attempt (expected $expected, got ${actual:-none})"
    rm -f "$gz"
  done
  [ "$actual" = "$expected" ] || { echo "ERROR: could not download $gz intact"; exit 1; }
  gunzip -t "$gz"
  gunzip -c "$gz" > "$FAMDB/$stem.part"
  mv "$FAMDB/$stem.part" "$FAMDB/$stem"
  printf '%s\t%s\t%s\t%s\n' "$gz" "$expected" "$(stat -c %s "$gz")" \
         "$(stat -c %s "$FAMDB/$stem")" > "$MD5DIR/$gz.md5ok"
  rm -f "$gz" "$gz.md5"
done

# Every listed file must have landed.
for gz in $FILES; do
  [ -s "$FAMDB/${gz%.gz}" ] || { echo "ERROR: missing after download: ${gz%.gz}"; exit 1; }
done

# Prove FamDB can open the set, and that taxonomy lookups work.
famdb.py -i "$FAMDB" info | tee "$DEST/famdb_info.txt"
grep -q "Version" "$DEST/famdb_info.txt" || { echo "ERROR: famdb.py info failed"; exit 1; }
famdb.py -i "$FAMDB" lineage -ad -c --format totals "Drosophila melanogaster" \
  | tee "$DEST/famdb_smoke.txt"
grep -q "entries in ancestors" "$DEST/famdb_smoke.txt" \
  || { echo "ERROR: famdb.py taxonomy smoke test failed"; exit 1; }

# Provenance manifest (shown in every report).
{
  printf 'file\tmd5_gz\tbytes_gz\tbytes_h5\n'
  cat "$MD5DIR"/*.md5ok
  echo "source: $BASE_URL"
  grep -E "Version|Date|Database" "$DEST/famdb_info.txt" | head -5
  echo "installed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$DEST/famdb_manifest.tsv"

cd /
rm -rf "$WORK"
echo "== reference data contents"
ls -la "$FAMDB"
du -sh "$FAMDB"
cat "$DEST/famdb_manifest.tsv"
echo "DFAM_DOWNLOAD_OK"
