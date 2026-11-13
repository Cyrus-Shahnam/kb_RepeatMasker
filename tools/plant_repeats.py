#!/usr/bin/env python3
"""
Plant known repeat copies into a genome, for known-truth testing of masking.

Copies of library sequences are written OVER background sequence (replacement,
not insertion), so coordinates and contig lengths are unchanged and the
unplanted genome is an exact control.  Every copy is mutated with its own
seeded RNG, so a copy is identical whatever else is planted.

Library use:   from plant_repeats import plant, read_fasta, write_fasta
CLI (real-data experiment, lessons §12):
  python3 tools/plant_repeats.py --genome MG1655.fna --library dmel_dfam.fa \
      --copies 40 --divergence 0.02,0.05,0.10,0.15 --seed 7 \
      --out-fasta MG1655_planted.fna --manifest MG1655_planted.json
"""
import argparse
import json
import random
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    """Return [(id, seq)] in file order (whole file in memory; fixtures only)."""
    recs, sid, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if sid is not None:
                    recs.append((sid, "".join(chunks)))
                sid, chunks = line[1:].split()[0] if line[1:].split() else "", []
            elif line.strip():
                chunks.append(line.strip())
    if sid is not None:
        recs.append((sid, "".join(chunks)))
    return recs


def write_fasta(path, recs, width=60):
    with open(path, "w") as fh:
        for sid, seq in recs:
            fh.write(">%s\n" % sid)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def revcomp(s):
    return s.translate(COMP)[::-1]


def random_seq(rng, n, gc=0.42):
    at, cg = (1 - gc) / 2, gc / 2
    return "".join(rng.choices("ACGT", weights=(at, cg, cg, at), k=n))


def mutate(seq, rng, divergence, indel_fraction=0.15):
    """Substitutions + small indels at total rate `divergence`."""
    out = []
    for base in seq:
        r = rng.random()
        if r < divergence * (1 - indel_fraction):
            out.append(rng.choice([b for b in "ACGT" if b != base.upper()] or ["A"]))
        elif r < divergence * (1 - indel_fraction / 2):
            continue  # deletion
        elif r < divergence:
            out.append(base)
            out.append(rng.choice("ACGT"))  # insertion
        else:
            out.append(base)
    return "".join(out)


def plant(genome, library, copies, divergences, seed, spacing=1000,
          fragment_min=0.5, full_length_fraction=0.6, simple_repeats=()):
    """
    genome: [(id, seq)]; library: [(name, seq)] (name may contain '#class').
    Returns (planted_genome, manifest_dict).  Coordinates are 1-based inclusive.
    simple_repeats: [(unit, total_len)] planted as kind 'simple'.
    """
    genome = [(cid, list(seq)) for cid, seq in genome]
    occupied = {cid: [] for cid, _ in genome}
    lengths = {cid: len(s) for cid, s in genome}
    place_rng = random.Random("%s:place" % seed)
    plants = []

    def find_slot(length):
        for _ in range(5000):
            weights = [max(0, lengths[c] - length - 2 * spacing) for c, _ in genome]
            if not any(weights):
                break
            cid = place_rng.choices([c for c, _ in genome], weights=weights)[0]
            start = place_rng.randint(spacing, lengths[cid] - length - spacing)
            end = start + length - 1
            if all(end + spacing < a or start > b + spacing for a, b in occupied[cid]):
                occupied[cid].append((start, end))
                return cid, start, end
        raise ValueError("No room to plant a %d bp copy; use fewer/shorter copies"
                         % length)

    def write(cid, start, seq):
        g = dict(genome)[cid]
        g[start - 1:start - 1 + len(seq)] = list(seq)

    for i in range(copies):
        rng = random.Random("%s:copy:%d" % (seed, i))
        name, cons = library[i % len(library)]
        # every family meets every divergence level (no family/divergence confounding)
        div = divergences[(i // len(library)) % len(divergences)]
        if rng.random() < full_length_fraction or len(cons) < 400:
            frag = cons
            fb, fe = 1, len(cons)
        else:
            flen = int(len(cons) * rng.uniform(fragment_min, 0.9))
            fb = rng.randint(1, len(cons) - flen + 1)
            frag = cons[fb - 1:fb - 1 + flen]
            fe = fb + flen - 1
        copy = mutate(frag.upper(), rng, div)
        strand = "+" if rng.random() < 0.5 else "-"
        if strand == "-":
            copy = revcomp(copy)
        cid, start, end = find_slot(len(copy))
        write(cid, start, copy)
        family, _, rclass = name.partition("#")
        plants.append({"id": "plant%03d" % (i + 1), "kind": "te", "family": family,
                       "class": rclass or "Unknown", "contig": cid, "start": start,
                       "end": end, "strand": strand, "divergence": div,
                       "consensus_begin": fb, "consensus_end": fe,
                       "length": len(copy)})
    for j, (unit, total) in enumerate(simple_repeats):
        seq = (unit * (total // len(unit) + 1))[:total]
        cid, start, end = find_slot(len(seq))
        write(cid, start, seq)
        plants.append({"id": "simple%02d" % (j + 1), "kind": "simple",
                       "family": "(%s)n" % unit, "class": "Simple_repeat",
                       "contig": cid, "start": start, "end": end, "strand": "+",
                       "divergence": 0.0, "length": len(seq)})
    out = [(cid, "".join(seq)) for cid, seq in genome]
    manifest = {"seed": seed, "copies": copies, "divergences": list(divergences),
                "contigs": [{"id": c, "length": lengths[c]} for c, _ in genome],
                "plants": sorted(plants, key=lambda p: (p["contig"], p["start"]))}
    return out, manifest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--library", required=True,
                    help="FASTA of repeat consensus sequences (name#class/family)")
    ap.add_argument("--copies", type=int, default=40)
    ap.add_argument("--divergence", default="0.02,0.05,0.10,0.15")
    ap.add_argument("--seed", default="20261020")
    ap.add_argument("--min-length", type=int, default=300)
    ap.add_argument("--max-length", type=int, default=6000)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--manifest", required=True)
    a = ap.parse_args(argv)
    lib = [(n, s) for n, s in read_fasta(a.library)
           if a.min_length <= len(s) <= a.max_length]
    if not lib:
        sys.exit("No library sequences between %d and %d bp" % (a.min_length,
                                                                a.max_length))
    random.Random("%s:lib" % a.seed).shuffle(lib)
    divs = [float(x) for x in a.divergence.split(",")]
    out, man = plant(read_fasta(a.genome), lib, a.copies, divs, a.seed)
    man.update({"genome": a.genome, "library": a.library})
    write_fasta(a.out_fasta, out)
    with open(a.manifest, "w") as fh:
        json.dump(man, fh, indent=1)
    print("planted %d copies (%d bp) into %s -> %s" % (
        len(man["plants"]), sum(p["length"] for p in man["plants"]), a.genome,
        a.out_fasta))


if __name__ == "__main__":
    main()
