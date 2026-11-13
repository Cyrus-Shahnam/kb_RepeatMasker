#!/usr/bin/env python3
"""
FAKE famdb.py for dry runs (see tools/fake_bin/RepeatMasker). Mimics the
output lines of FamDB 3.0.0 that kb_repeatmasker parses. NEVER used in the image.
"""
import os
import sys

KNOWN = ("drosophila melanogaster", "fakeus testus")


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: famdb.py [-h] [-l LOG_LEVEL] [-i DB_DIR]\nThis is famdb.py "
              "version 3.0.0.")
        return 0
    args = list(argv)
    if args[0] == "-i":
        args = args[2:]
    cmd, rest = args[0], args[1:]
    term = rest[-1] if rest else ""
    if cmd == "info":
        print("FamDB Directory               : fake\nFamDB Creation Format Version :"
              " 3.0.0\n\nDatabase : Dfam\nVersion  : 4.0\nDate     : 2026-06-01\n")
        return 0
    if term.lower() not in KNOWN:
        sys.stderr.write("No species found for search term '%s'\n" % term)
        return 0  # the real famdb.py also exits 0 here
    if cmd == "lineage":
        print("120 entries in ancestors; 30 lineage-specific entries; 120 ancestral "
              "entries present; 30 lineage-specific entries present")
    elif cmd == "check":
        print("Partition check for '%s' (tax id: 7227):\n" % term)
        print("  Curated Consensus    partition 0 [root]:   present")
        if os.environ.get("FAKE_FAMDB_MISSING"):
            print("  Uncurated Consensus  partition 9 [Insecta]: MISSING  "
                  "[dfam40.uncurated.consensus.9.h5]")
    elif cmd == "families":
        lib = os.environ.get("FAKE_SPECIES_LIB")
        if lib:
            for line in open(lib):
                if line.startswith(">"):
                    line = line.rstrip("\n") + " @Drosophila [S:40]\n"
                sys.stdout.write(line)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
