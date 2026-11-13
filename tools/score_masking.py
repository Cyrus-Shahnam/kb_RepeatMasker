#!/usr/bin/env python3
"""
Score a masked FASTA against a plant_repeats.py manifest.

Masked = lowercase, N or X in the masked FASTA (works for soft and hard masking).
Pass --unmasked (the planted FASTA before masking) when the input had its own N
runs or lowercase, so pre-existing masking is not counted.

  python3 tools/score_masking.py --manifest planted.json --masked planted.masked.fa \
      [--unmasked planted.fna] [--control-masked control.masked.fa]
"""
import argparse
import json
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from plant_repeats import read_fasta  # noqa: E402


def mask_vector(seq, orig=None):
    if orig is None:
        return [c.islower() or c in "NnXx" for c in seq]
    return [(c.islower() and not o.islower()) or (c in "NnXx" and o not in "NnXx")
            for c, o in zip(seq, orig)]


def score(manifest, masked_recs, unmasked_recs=None, detect_fraction=0.5):
    masked = dict(masked_recs)
    orig = dict(unmasked_recs or [])
    vec = {}
    for c in manifest["contigs"]:
        if c["id"] not in masked:
            raise ValueError("contig %s missing from masked FASTA" % c["id"])
        if len(masked[c["id"]]) != c["length"]:
            raise ValueError("contig %s length %d != manifest %d"
                             % (c["id"], len(masked[c["id"]]), c["length"]))
        vec[c["id"]] = mask_vector(masked[c["id"]], orig.get(c["id"]))
    per_plant, inside = [], {c: set() for c in vec}
    for p in manifest["plants"]:
        v = vec[p["contig"]][p["start"] - 1:p["end"]]
        frac = sum(v) / float(len(v))
        inside[p["contig"]].update(range(p["start"] - 1, p["end"]))
        per_plant.append(dict(p, masked_fraction=round(frac, 4),
                              detected=frac >= detect_fraction))
    te = [p for p in per_plant if p["kind"] == "te"]
    te_bp = sum(p["length"] for p in te)
    te_masked = sum(p["masked_fraction"] * p["length"] for p in te)
    bg_total = bg_masked = 0
    for cid, v in vec.items():
        for i, x in enumerate(v):
            if i not in inside[cid]:
                bg_total += 1
                bg_masked += x
    by_div = {}
    for p in te:
        d = by_div.setdefault(str(p["divergence"]), [0, 0])
        d[0] += p["detected"]
        d[1] += 1
    return {
        "te_copies": len(te), "te_detected": sum(p["detected"] for p in te),
        "te_bp": te_bp, "te_bp_recall": round(te_masked / max(1, te_bp), 4),
        "simple": [p for p in per_plant if p["kind"] == "simple"],
        "background_bp": bg_total, "background_masked_bp": bg_masked,
        "background_masked_frac": round(bg_masked / float(max(1, bg_total)), 6),
        "detected_by_divergence": {k: "%d/%d" % tuple(v) for k, v in sorted(by_div.items())},
        "per_plant": per_plant,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--masked", required=True)
    ap.add_argument("--unmasked")
    ap.add_argument("--control-masked",
                    help="masked FASTA of the UNPLANTED genome (the control)")
    ap.add_argument("--json-out")
    a = ap.parse_args(argv)
    man = json.load(open(a.manifest))
    masked_recs = read_fasta(a.masked)
    res = score(man, masked_recs, read_fasta(a.unmasked) if a.unmasked else None)
    if a.control_masked:
        ctrl = {c: mask_vector(s) for c, s in read_fasta(a.control_masked)}
        masked = dict(masked_recs)
        res["control_masked_bp"] = sum(sum(v) for v in ctrl.values())
        # planted loci that were ALREADY masked in the control are not attributable
        # to the planted copy: report recall over control-unmasked bp as well.
        new = tot = 0
        for p in res["per_plant"]:
            if p["kind"] != "te":
                continue
            cv = ctrl[p["contig"]][p["start"] - 1:p["end"]]
            mv = mask_vector(masked[p["contig"]][p["start"] - 1:p["end"]])
            for c_, m_ in zip(cv, mv):
                if not c_:
                    tot += 1
                    new += m_
        res["te_bp_recall_vs_control"] = round(new / float(max(1, tot)), 4)
    print("TE copies detected: %d/%d   bp recall: %.3f   background masked: %d bp "
          "(%.4f%%)" % (res["te_detected"], res["te_copies"], res["te_bp_recall"],
                        res["background_masked_bp"],
                        100 * res["background_masked_frac"]))
    print("detected by divergence:", res["detected_by_divergence"])
    for s in res["simple"]:
        print("simple %s masked fraction %.2f" % (s["family"], s["masked_fraction"]))
    if "control_masked_bp" in res:
        print("control masked bp: %d   recall vs control: %.3f"
              % (res["control_masked_bp"], res["te_bp_recall_vs_control"]))
    if a.json_out:
        with open(a.json_out, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
