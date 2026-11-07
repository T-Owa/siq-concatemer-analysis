# Reviewer guide: SIQ concatemer junction analysis

Software version: **0.1.0**. Repository: [T-Owa/siq-concatemer-analysis](https://github.com/T-Owa/siq-concatemer-analysis).

This archive supports numerical reproduction of the supplied **SIQ event-table post-processing**. It includes actual processed study data, the reference configuration, sample metadata and independent historical comparison tables. Commands below are run from the repository root.

## Inputs and scope

| Dataset | Samples | SIQ event rows | Sum of reported event counts |
|---|---:|---:|---:|
| Donor modifications (`modify`) | 60 | 43,532 | 1,992,656 |
| Repair inhibitors (`inhibitor`) | 54 | 27,771 | 1,997,053 |

Each condition/time point has three independent biological experiments. The analysis reference is a **196-nt head-to-tail amplicon from the 825-bp donor**, with junction index 136 in zero-based coordinates. Reference sequence and primers are in `config/junction.json`.

The gzip CSV inputs retain all SIQ-reported rows, including zero-count rows. These are processed SIQ event tables, not FASTQ reads. This archive does not rerun upstream read merging or SIQ and cannot recover events excluded by those upstream steps. See [Data and outputs](DATA.md) for fields and provenance.

## Recommended numerical reproduction

Python **3.10 or newer** is sufficient for the following commands; no third-party Python packages are required. On Windows, use `python` instead of `python3` where appropriate.

```bash
python3 scripts/run_analysis.py --outdir review-results
python3 scripts/verify_reference_results.py --results review-results --report review-results/verification.json
python3 -m unittest discover -s tests -v
```

Use a new or empty output directory. The workflow runs both datasets by default. Successful verification reports `status: passed`, both dataset names, and **24,531 checked values**. All **32 unit tests** should pass.

Reference comparison uses exact equality for count/identifier checks and `math.isclose` with relative and absolute tolerances of `1e-10` for floating-point checks. This verifies agreement with the archived study analyses, not an independent biological validation of the experiment.

## Analysis-to-output map

Paths below are relative to `review-results/<dataset>/`, except the run manifest.

| Analysis | Output to inspect | Population / denominator |
|---|---|---|
| Event composition | `summary/event_types_by_sample.csv`, `summary/event_types_by_group.csv` | All SIQ-reported `countEvents` for each sample |
| Deletion length and position | `summary/deletion_lengths_by_sample.csv`, `summary/deletion_lengths_by_group.csv`, `summary/deletion_positions_by_sample.csv`, `summary/deletion_positions_by_group.csv` | Pure `Type == DELETION`; alternative all-event fractions are also supplied |
| Strict junction MH | `mh/events.csv`, `summary/mh_distribution_by_sample.csv`, `summary/mh_status_by_sample.csv` | Pure deletions with exact full-Raw reconstruction and a primer-preserving junction-compatible placement |
| All-length dominant patterns | `patterns/pattern_definitions.csv`, `patterns/all_top5.csv` | All pure deletion counts in each sample |
| Modification-only 6/12-bp patterns | `patterns/length_sample_totals.csv`, `patterns/length_group_patterns.csv`, `patterns/length_top5.csv`, `patterns/length_coverage.csv` | All pure deletions of the specified length, retaining internal-only and Raw-disagreement records |
| Matched-length/descriptive MH comparisons | `summary/exact_length_mh_by_sample.csv`, `summary/exact_length_mh_by_group.csv`, `summary/length_standardized.csv` (modifications) | Strict-MH-evaluable deletions at the stated lengths; common-length weighting as defined in Methods |

`summary/sample_metrics.csv` and `summary/group_metrics.csv` collect metrics from several analyses; each metric retains its analysis-specific population and denominator.

`review-results/run_manifest.json` records input/configuration hashes, Python version and executed analysis steps. Historical checks are in `expected/`; their source identities are described in `expected/fixture_provenance.json`.

## Optional figure inspection

After numerical reproduction, install the optional dependencies using the [README](../README.md#quick-start) instructions, then run:

```bash
python3 scripts/plot_figures.py --results review-results --outdir review-results/figures
```

The default plotting command produces five modification-dataset figures for exact-length MH and 6/12-bp deletion patterns, as editable SVG, PDF and PNG, with captions and plotting metadata. It does not render every analysis table or assign manuscript figure numbers. Plotting does not change the calculated values.

## Interpretation and verification scope

- Fractions are calculated per sample, then averaged with equal weight across biological experiments; SD is sample SD. Reads are not biological replicates.
- Strict MH excludes 6,061 modification and 1,428 inhibitor deletion reads. Unevaluable MH is NA, not zero.
- The 6-bp and 12-bp pattern denominators contain 78,572 and 47,634 reads. They differ from strict-MH denominators and are not renormalized to the displayed Top5.
- MH is a maximum sequence-compatible candidate; it does not demonstrate use of MMEJ/TMEJ or another repair pathway. Equivalent coordinates do not uniquely identify cleavage positions.
- No inferential statistical tests are performed. The length-standardized comparisons are descriptive.

Full definitions are in [Methods](METHODS.md). [Validation](VALIDATION.md) separates the original validation record from the submission-package rerun. `docs/SUBMISSION_VALIDATION.json` records the latter, including the source and input identities.

On Linux/WSL, verify the extracted archive contents before generating outputs:

```bash
sha256sum -c SHA256SUMS
```
