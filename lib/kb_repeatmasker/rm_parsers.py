"""
Pure-Python helpers for kb_repeatmasker: FASTA streaming, contig-ID mapping,
RepeatMasker .out / .tbl / .gff parsing and masking statistics.

Nothing in this file talks to KBase services, so it can be unit-tested (and
dry-run) outside the container.
"""
import re

# RepeatMasker's FastaDB truncates/rejects IDs longer than 50 characters, and it
# refuses file names with spaces or shell characters.  We therefore always rename
# contigs to short safe IDs before running and translate the outputs back.
SHORT_ID_FMT = "rmq%06d"

# Repeat classes that are NOT interspersed (transposable-element derived) repeats.
NON_INTERSPERSED = ("Simple_repeat", "Low_complexity", "Satellite", "Small RNA")
CATEGORY_ORDER = ("SINE", "LINE", "LTR", "DNA", "RC", "Retroposon", "Unknown",
                  "Other", "Satellite", "Small RNA", "Simple_repeat",
                  "Low_complexity")


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #
def iter_fasta(path):
    """Yield (seq_id, description, sequence) one record at a time."""
    seq_id, desc, chunks = None, "", []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, desc, "".join(chunks)
                header = line[1:].strip()
                parts = header.split(None, 1)
                seq_id = parts[0] if parts else ""
                desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            elif line.strip():
                if seq_id is None:
                    raise ValueError("%s: sequence data before the first '>' header "
                                     "(is the file gzipped or not FASTA?)" % path)
                chunks.append(line.strip())
    if seq_id is not None:
        yield seq_id, desc, "".join(chunks)


def write_fasta_record(fh, seq_id, seq, width=60):
    fh.write(">%s\n" % seq_id)
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def rename_fasta(in_path, out_path):
    """
    Copy a FASTA with contigs renamed to short IDs.
    Returns a list of dicts: short_id, original_id, length (input order).
    """
    mapping = []
    seen = set()
    with open(out_path, "w") as out:
        for n, (sid, _desc, seq) in enumerate(iter_fasta(in_path), start=1):
            if not sid:
                raise ValueError("Empty FASTA header in %s (record %d)" % (in_path, n))
            if sid in seen:
                raise ValueError("Duplicate contig ID %r in %s" % (sid, in_path))
            seen.add(sid)
            short = SHORT_ID_FMT % n
            write_fasta_record(out, short, seq)
            mapping.append({"short_id": short, "original_id": sid, "length": len(seq)})
    if not mapping:
        raise ValueError("No sequences found in %s" % in_path)
    return mapping


def write_id_map(mapping, path):
    with open(path, "w") as fh:
        fh.write("short_id\toriginal_id\tlength\n")
        for m in mapping:
            fh.write("%s\t%s\t%d\n" % (m["short_id"], m["original_id"], m["length"]))


def gc_and_length(path):
    total = gc = acgt = 0
    for _sid, _d, seq in iter_fasta(path):
        s = seq.upper()
        total += len(s)
        g = s.count("G") + s.count("C")
        gc += g
        acgt += g + s.count("A") + s.count("T")
    return total, (100.0 * gc / acgt if acgt else 0.0)


# --------------------------------------------------------------------------- #
# Masked FASTA: translate IDs back and measure masking
# --------------------------------------------------------------------------- #
def count_masked(orig_seq, masked_seq, mask_type):
    """Number of positions RepeatMasker masked in one contig."""
    if mask_type == "soft":
        return sum(1 for c in masked_seq if c.islower())
    if mask_type == "hard_X":
        return masked_seq.upper().count("X") - orig_seq.upper().count("X")
    # hard_N: RepeatMasker only ever turns bases into N, never back.
    return masked_seq.upper().count("N") - orig_seq.upper().count("N")


def translate_masked_fasta(input_fasta, masked_fasta, out_fasta, mapping, mask_type):
    """
    Walk the original and RepeatMasker-masked FASTA in step (RepeatMasker keeps
    input order), write the masked sequences under their ORIGINAL contig IDs,
    and return per-contig masking counts plus integrity problems found.
    """
    problems = []
    per_contig = []
    by_short = {m["short_id"]: m for m in mapping}
    orig_iter = iter_fasta(input_fasta)
    seen = set()
    with open(out_fasta, "w") as out:
        for short, _d, mseq in iter_fasta(masked_fasta):
            try:
                oid, _od, oseq = next(orig_iter)
            except StopIteration:
                problems.append("masked FASTA has more records than the input")
                break
            m = by_short.get(short)
            if m is None:
                problems.append("unexpected record %r in masked FASTA" % short)
                continue
            if m["original_id"] != oid:
                problems.append("record order changed at %s (expected %s)"
                                % (short, oid))
            if len(mseq) != len(oseq):
                problems.append("length changed for %s: %d -> %d"
                                % (oid, len(oseq), len(mseq)))
            seen.add(short)
            n_masked = count_masked(oseq, mseq, mask_type)
            per_contig.append({"contig_id": m["original_id"], "length": len(mseq),
                               "masked_bp": n_masked})
            write_fasta_record(out, m["original_id"], mseq)
    missing = [m["original_id"] for m in mapping if m["short_id"] not in seen]
    if missing:
        problems.append("%d contig(s) missing from masked output (first: %s)"
                        % (len(missing), missing[0]))
    return per_contig, problems


def copy_unmasked(input_fasta, out_fasta):
    """Used when RepeatMasker found nothing and wrote no .masked file."""
    per_contig = []
    with open(out_fasta, "w") as out:
        for sid, _d, seq in iter_fasta(input_fasta):
            write_fasta_record(out, sid, seq)
            per_contig.append({"contig_id": sid, "length": len(seq), "masked_bp": 0})
    return per_contig


# --------------------------------------------------------------------------- #
# .out (the primary per-annotation table)
# --------------------------------------------------------------------------- #
_OUT_ROW = re.compile(r"^(\s*\d+\s+\S+\s+\S+\s+\S+\s+)(\S+)(.*)$")
NO_REPEATS_MARKERS = ("There were no repetitive sequences detected",
                      "No repetitive sequences were detected")


def parse_out(path, id_map=None):
    """
    Parse a RepeatMasker .out file.  Returns (rows, no_repeats_flag).
    Columns: score div del ins query qbegin qend (left) strand repeat class/family
             rbegin/(left) rend rleft/begin ID [*]
    """
    id_map = id_map or {}
    rows, no_repeats = [], False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if any(mk in s for mk in NO_REPEATS_MARKERS):
                no_repeats = True
                continue
            # "(  123)" -> "(123)": a padded "(left)" column must stay ONE token,
            # otherwise every later column shifts (strand read as repeat name).
            f = re.sub(r"\(\s+", "(", s).split()
            if not f[0].isdigit() or len(f) < 14:
                continue  # header lines ("SW", "score") and anything unexpected
            try:
                qbeg, qend = int(f[5]), int(f[6])
            except ValueError:
                continue
            rows.append({
                "score": int(f[0]), "div": float(f[1]), "del": float(f[2]),
                "ins": float(f[3]),
                "contig": id_map.get(f[4], f[4]),
                "begin": min(qbeg, qend), "end": max(qbeg, qend),
                "strand": "-" if f[8] == "C" else "+",
                "repeat": f[9], "class": f[10],
                "id": f[14] if len(f) > 14 else "",
                "overlapped": f[-1] == "*",
            })
    return rows, no_repeats


def translate_out_file(in_path, out_path, id_map):
    """Rewrite the query-sequence column of a .out file with original IDs."""
    with open(in_path) as src, open(out_path, "w") as dst:
        for line in src:
            m = _OUT_ROW.match(line.rstrip("\n"))
            if m and m.group(2) in id_map:
                dst.write(m.group(1) + id_map[m.group(2)] + m.group(3) + "\n")
            else:
                dst.write(line)


def translate_gff_file(in_path, out_path, id_map):
    with open(in_path) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith("#"):
                toks = line.rstrip("\n").split(" ")
                dst.write(" ".join(id_map.get(t, t) for t in toks) + "\n")
                continue
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] in id_map:
                cols[0] = id_map[cols[0]]
            dst.write("\t".join(cols) + "\n")


def repeat_category(rclass):
    top = rclass.split("/")[0]
    if top in ("Simple_repeat", "Low_complexity", "Satellite", "SINE", "LINE",
               "LTR", "DNA", "RC", "Retroposon"):
        return top
    if top == "PLE":
        return "LINE"
    if "RNA" in top:
        return "Small RNA"
    if top in ("Unknown", "Unspecified", ""):
        return "Unknown"
    return "Other"


def summarize_rows(rows):
    """Aggregate annotations per category and per family (annotated bp)."""
    cats, fams = {}, {}
    for r in rows:
        length = r["end"] - r["begin"] + 1
        c = repeat_category(r["class"])
        e = cats.setdefault(c, {"count": 0, "bp": 0})
        e["count"] += 1
        e["bp"] += length
        key = (r["repeat"], r["class"])
        fe = fams.setdefault(key, {"repeat": r["repeat"], "class": r["class"],
                                   "count": 0, "bp": 0})
        fe["count"] += 1
        fe["bp"] += length
    interspersed = sum(v["bp"] for k, v in cats.items() if k not in NON_INTERSPERSED)
    simple_low = sum(cats.get(k, {"bp": 0})["bp"]
                     for k in ("Simple_repeat", "Low_complexity"))
    top = sorted(fams.values(), key=lambda x: (-x["bp"], x["repeat"]))
    return {"categories": cats, "interspersed_bp": interspersed,
            "simple_lowcomp_bp": simple_low, "top_families": top[:15],
            "n_families": len(fams)}


# --------------------------------------------------------------------------- #
# .tbl (summary; used only as a cross-check)
# --------------------------------------------------------------------------- #
def parse_tbl(path):
    with open(path) as fh:
        text = fh.read()
    out = {}
    m = re.search(r"total length:\s+(\d+) bp", text)
    if m:
        out["total_length"] = int(m.group(1))
    m = re.search(r"bases masked:\s+(\d+) bp\s*\(\s*([\d.]+)\s*%\)", text)
    if m:
        out["bases_masked"] = int(m.group(1))
        out["pct_masked"] = float(m.group(2))
    m = re.search(r"GC level:\s+([\d.]+)\s*%", text)
    if m:
        out["gc_level"] = float(m.group(1))
    m = re.search(r"sequences:\s+(\d+)", text)
    if m:
        out["sequences"] = int(m.group(1))
    return out


# --------------------------------------------------------------------------- #
# famdb.py output
# --------------------------------------------------------------------------- #
_TOTALS = re.compile(r"(\d+) entries in ancestors; (\d+) lineage-specific entries"
                     r"(?:; (\d+) ancestral entries present; (\d+) lineage-specific"
                     r" entries present)?")


def parse_lineage_totals(text):
    m = _TOTALS.search(text)
    if not m:
        return None
    vals = [int(v) if v is not None else None for v in m.groups()]
    return {"ancestral": vals[0], "lineage_specific": vals[1],
            "ancestral_present": vals[2], "lineage_specific_present": vals[3]}


def parse_famdb_check(text):
    """Return the list of MISSING partition file names from `famdb.py check`."""
    missing = []
    for line in text.splitlines():
        if "MISSING" in line:
            # "partition 50 [Eutheria]:  MISSING  [dfam40.uncurated.hmm.50.h5]"
            # - the file name is the bracket AFTER "MISSING", not the taxon one.
            m = re.search(r"MISSING\s*\[([^\]]+)\]", line)
            missing.append(m.group(1) if m else line.strip())
    return missing


def parse_famdb_info(text):
    info = {}
    for key, pat in (("db_version", r"^\s*Version\s*:\s*(\S+)"),
                     ("db_date", r"^\s*Date\s*:\s*(\S+)"),
                     ("db_name", r"^\s*Database\s*:\s*(\S+)")):
        m = re.search(pat, text, re.M)
        if m:
            info[key] = m.group(1)
    return info


def parse_rm_version(text):
    m = re.search(r"RepeatMasker version\s+(\S+)", text)
    return m.group(1) if m else None


def count_fasta_records(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n
