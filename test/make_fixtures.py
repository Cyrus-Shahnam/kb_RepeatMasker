#!/usr/bin/env python3
"""
Seeded known-truth fixtures for kb_repeatmasker (lessons §8).

  python3 test/make_fixtures.py <out_dir>

Writes:
  te_library.fa      4 synthetic TE consensus sequences, IDs "name#class/family"
  genome_control.fa  random background, 3 contigs (one ID > 50 characters)
  genome_planted.fa  the same background with 16 diverged TE copies and
                     2 microsatellites written OVER it (same coordinates)
  manifest.json      exactly what was planted where

Nothing is committed except this script; tests regenerate fixtures in scratch.
"""
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from plant_repeats import plant, random_seq, write_fasta  # noqa: E402

SEED = "20261020"
LONG_ID = ("contig_with_a_deliberately_long_identifier_from_an_assembler_"
           "NODE_1_length_70000_cov_31.2")  # > 50 chars: RepeatMasker limit
CONTIGS = [(LONG_ID, 70000), ("ctg2", 50000), ("ctg3", 30000)]
TE_FAMILIES = [("kbTE_LINE1#LINE/L1", 3000), ("kbTE_LTR1#LTR/Gypsy", 1500),
               ("kbTE_DNA1#DNA/hAT", 900), ("kbTE_SINE1#SINE/tRNA", 320)]
DIVERGENCES = (0.02, 0.05, 0.08, 0.12)
N_COPIES = 16  # 4 families x 4 divergence levels
SIMPLE = [("CA", 120), ("TTAGGG", 96)]


def build(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    background = [(cid, random_seq(random.Random("%s:bg:%s" % (SEED, cid)), n))
                  for cid, n in CONTIGS]
    library = [(name, random_seq(random.Random("%s:te:%s" % (SEED, name)), n, gc=0.40))
               for name, n in TE_FAMILIES]
    planted, manifest = plant(background, library, N_COPIES, DIVERGENCES, SEED,
                              simple_repeats=SIMPLE)
    manifest["fixture"] = "kb_repeatmasker synthetic known-truth"
    paths = {k: os.path.join(out_dir, v) for k, v in (
        ("library", "te_library.fa"), ("control", "genome_control.fa"),
        ("planted", "genome_planted.fa"), ("manifest", "manifest.json"))}
    write_fasta(paths["library"], library)
    write_fasta(paths["control"], background)
    write_fasta(paths["planted"], planted)
    with open(paths["manifest"], "w") as fh:
        json.dump(manifest, fh, indent=1)
    return paths, manifest


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data", "fixtures")
    p, man = build(out)
    te = [x for x in man["plants"] if x["kind"] == "te"]
    print("wrote %s: %d TE copies (%d bp), %d simple repeats, contigs %s" % (
        out, len(te), sum(x["length"] for x in te), len(man["plants"]) - len(te),
        ", ".join("%s:%d" % (c["id"][:12], c["length"]) for c in man["contigs"])))
