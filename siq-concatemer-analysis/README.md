# SIQ concatemer junction analysis

Reproducible post-processing of SIQ-reported head-to-tail donor junctions: event composition, deletion length and position, sequence-compatible microhomology (MH), and dominant deletion patterns.

This repository contains the custom analysis used for an 825-bp linear dsDNA donor study. The supplied data comprise **60 modification samples and 54 inhibitor samples**, with three independent biological experiments per condition and time point (12 h and 24 h). It includes actual SIQ event tables, not simulated counts. SIQ itself and the sequencing FASTQ files are not included.

**MH is a maximum sequence-compatible candidate, not evidence that a particular repair pathway was used.** See [the analysis definitions](docs/METHODS.md) before interpreting the outputs.

Japanese instructions: [README_JA.md](README_JA.md).

## Quick start

Python 3.10 or newer is required. The bundled CSV analysis uses the Python standard library. Run the commands from this repository's top-level folder:

```bash
python3 scripts/run_analysis.py --outdir results
python3 scripts/verify_reference_results.py --results results
python3 -m unittest discover -s tests -v
```

Choose a new or empty output directory. Outputs are not overwritten. On Windows PowerShell, use `python` in place of `python3` when appropriate. Under Ubuntu/WSL, use `python3`, not `py`.

To generate editable SVG, PDF and PNG figures, install the optional dependencies and use a new output directory. The pinned dependencies require Python 3.11 or newer; this package was validated with Python 3.12.14:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 scripts/run_analysis.py --outdir results_with_figures --figures
```

Alternatively, add figures to an existing analysis:

```bash
python3 scripts/plot_figures.py --results results --outdir results/figures
```

These plots are a portable static plotting implementation. They do not recreate the Figure Studio web interface or guarantee pixel-identical formatting of earlier figures.

## Contents

| Location | Purpose |
|---|---|
| `scripts/analyze_junction_mh.py` | Preserved original MH algorithm and read-only SIQ Excel entry point |
| `scripts/import_siq.py` | Extract analysis fields from an SIQ workbook, checking RunInfo |
| `scripts/annotate_events.py` | Apply the same MH algorithm to portable CSV or gzip CSV inputs |
| `scripts/summarize_events.py` | Event, deletion length/position, strict MH and descriptive length-standardized summaries |
| `scripts/aggregate_patterns.py` | Allele grouping, all-length Top5 and modification-only 6/12-bp Top5 |
| `scripts/plot_figures.py` | Static figures from the calculated tables |
| `scripts/run_analysis.py` | Run the bundled data through the complete post-processing workflow |
| `scripts/verify_reference_results.py` | Compare regenerated values with archived analysis results |
| `config/` | Reference sequence, primers, junction and historical pattern ID mapping |
| `data/` | Selected rawData fields, sample manifest and source provenance |
| `expected/` | Archived numerical reference results and independent validation evidence |
| `tests/` | Scientific boundary and aggregation tests |

## Populations and denominators

All counts are weighted by `countEvents`. “All events” means the sum of counts in SIQ `rawData`, **not** MergedCorrect or all sequenced reads. SIQ output filters have already acted on these data.

| Analysis | Population / denominator |
|---|---|
| Event composition | All reported event counts for each sample |
| Deletion length and original-coordinate position coverage | Pure `Type == DELETION`; alternative denominator of all reported events is also supplied |
| Strict MH | Pure deletions with valid coordinates, exact full Raw reconstruction, and at least one primer-preserving junction-compatible placement |
| All-length Top5 | All pure deletion reads; coordinate-reconstructed alleles are grouped, with Raw disagreements flagged |
| 6/12-bp Top5 | All pure deletion reads of the specified length, including internal-only and Raw-disagreement records; modifications only |

Sample fractions are calculated first, then averaged with equal weight across independent experiments. SD is the sample SD. Read pooling does not replace biological replication. Missing denominator values are not converted to zero. Pattern absence in a sample with a positive denominator is counted as zero.

Within the supplied data, 6-bp and 12-bp Top5 denominators contain 78,572 and 47,634 reads, respectively. The strict MH analysis excludes 6,061 modification and 1,428 inhibitor deletion reads. The two analyses therefore must not share denominators silently.

## Other SIQ exports

The original Excel entry point can be used directly:

```bash
python3 scripts/analyze_junction_mh.py --input your_SIQ.xlsx --config config/junction.json --outdir mh_excel_results
```

It requires `openpyxl`. To export a portable input table:

```bash
python3 scripts/import_siq.py --input your_SIQ.xlsx --dataset modify --config config/junction.json --outdir imported
python3 scripts/run_analysis.py --dataset modify --data-dir imported --outdir new_results
```

The configuration is study-specific. Change it for another reference or junction. The sample manifest must correctly identify independent experiments. Pattern catalog IDs are specific to the supplied reference and must not be reused for unrelated sequences. See [data documentation](docs/DATA.md).

## Scope and provenance

Version 0.1.0 is a publication-preparation package assembled on 2026-11-04. The original MH algorithm is retained; portable import, aggregation, plotting and orchestration were consolidated for this package. Numerical reproduction is checked against the earlier archived analyses. The final manuscript figures and analysis selection remain to be linked to a release.

No inferential statistical tests are run in this package. The length-standardized MH comparisons are descriptive. This package does not reproduce qPCR concatemer formation rates, RAISING analyses, read merging, or the upstream SIQ computation.

SIQ 2.1 and the actual run settings were recovered from the source workbooks. Cite SIQ separately: [official SIQ repository](https://github.com/RobinVanSchendel/SIQ), [SIQ paper](https://doi.org/10.1093/nargab/lqac063). This project is a custom downstream analysis and is not an official SIQ distribution.

See [publication preparation](docs/PUBLICATION.md) for release metadata, licensing scope and code availability wording.
