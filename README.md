# kb_repeatmasker

KBase SDK module wrapping **RepeatMasker 4.2.4** (RMBlast engine + TRF) with the
**Dfam 4.0** curated library, installed as KBase reference data.

One app, `run_repeatmasker` ("Mask Repeats with RepeatMasker"):
Assembly (one or more) [+ optional custom repeat library Assembly] → masked
Assembly per input + report.

## What it does

| Library mode | When | RepeatMasker call |
|---|---|---|
| Dfam | species/taxon given | `-species "<taxon>"` (staged search, Dfam curated families) |
| Custom | only a library Assembly given | `-lib library.fa` |
| Combined | both given | `famdb.py families -ad --curated <taxon>` + library → `-lib` |

The taxon is resolved with `famdb.py lineage` **before** RepeatMasker runs, so a
misspelled name fails in seconds with a readable message, and the report states
how many Dfam families applied. `famdb.py check` warns when a needed partition is
not installed.

Outputs per input: masked Assembly (`<name>_masked`; soft = lowercase by default,
or hard N / X), and a zip with `<name>.out` (annotation table), `<name>.gff`,
`<name>.tbl`, `<name>.masked.fa`, `contig_id_map.tsv`, `per_contig_masking.tsv`
and the log. The report has repeat content by class, top families, warnings,
notes, parameters, exact commands, Dfam provenance and citations.

## Design decisions (and why)

- **Contig IDs are renamed** to `rmq000001…` before RepeatMasker and translated
  back in every output: RepeatMasker's FastaDB limits IDs to 50 characters and
  rejects file names with spaces/special characters. Tests use a >50-char ID.
- **Isolated conda env** (`/opt/conda/envs/repeatmasker`, micromamba, bioconda
  `repeatmasker=4.2.4`, `famdb=3.0.0`), never on the image PATH; injected into the
  PATH of each subprocess only.
- **FamDB data at `/data/famdb`** (read-only at runtime). The Dockerfile points
  FamDB there via `famdb.conf` and `Libraries/famdb → /data/famdb`.
  RepeatMasker calls `famdb.py info` on every run, so this data is needed even
  for custom-library-only runs.
- **Parameters that cannot take effect are dropped with a warning** (`cutoff`
  without a custom library, `uncurated` without a taxon, …).
- **No run-time partition download.** Partition sizes vary widely; add files to
  `data/refdata_famdb.txt` instead (the report names the exact file).
- `-pa` = CPUs / 4 (RMBlast uses 4 threads per batch). `-excln` is always on so
  scaffold N-runs do not dilute percentages.

## Reference data

`data/refdata_famdb.txt` lists the files (pinned to `Dfam_4.0`, not `current`):
the root partition `dfam40.0.h5` and `dfam40.curated.consensus.0.h5`.
`scripts/download_dfam_data.sh` (run by `entrypoint.sh init`) downloads each,
checks Dfam's `.md5` sidecar, decompresses, runs `famdb.py info` and a taxonomy
smoke test, and writes `/data/famdb_manifest.tsv`. After changing the list,
**bump `data-version` in `kbase.yml` and re-register**.

## Development

```bash
python3 tools/dry_run.py        # whole Runner with fake tools, no Docker (seconds)
make && kb-sdk validate          # in that order
python3 tools/check_names.py     # name chain + placeholders
kb-sdk test                      # builds image, runs init (downloads Dfam), tests
bash tools/diagnose.sh           # read-only state dump when something is wrong
```

Tests: fast tier (validation, Narrative-style blank/string params, command
order, parsers, fixture sanity) and a known-truth tier: 16 synthetic TE copies
(4 families × 2/5/8/12 % divergence) and 2 microsatellites planted over a seeded
random background, with the unplanted background as the control, plus a Dfam
tier (6 curated *D. melanogaster* families planted at 3 %). The Dfam tier is
**skipped** if `/data/famdb` is missing — "N passed, M skipped" means it did not run.

Real-data known truth: `tools/plant_repeats.py` plants Dfam consensus copies into
a real genome (control = the unplanted genome); `tools/score_masking.py` scores
the masked assemblies.

## Licences

RepeatMasker: Open Software License v2.1 (ISB). Dfam: CC0 1.0. TRF, RMBlast and
HMMER under their own licences (installed from bioconda). Cite RepeatMasker,
Dfam and TRF (see the report).
