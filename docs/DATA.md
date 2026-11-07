# Data and outputs

## Portable source tables

`data/modify_events.csv.gz` and `data/inhibitor_events.csv.gz` contain selected fields from the original SIQ `rawData` sheets. They retain every reported row, including zero-count reference rows. gzip CSV files can be read directly by the supplied scripts. Decompression is not required.

| Field | Meaning |
|---|---|
| source_excel_row | Original 1-based Excel row number |
| Alias | SIQ sample name, matched to the sample manifest |
| Type | Original SIQ event type |
| countEvents | Reported read count, used as a weight |
| Raw | Entire SIQ-reported event sequence |
| delStart, delEnd, delSize | Original SIQ deletion coordinates and length |
| homology, homologyLength | Original SIQ annotations, retained for provenance; not the custom MH result |

`data/samples.csv` identifies dataset, Alias, condition, condition_label, time and replicate. Source aliases are preserved. `Bractppin` is displayed as Bractoppin and `DI03` as D-I03. Modification aliases ds, phos, Bio, dsp and sp18 are displayed as Unmodified, Phos, Biotin, dSp and Sp18.

| Dataset | Source rows | Samples | Reported reads | Pure deletion reads | Strict-MH-evaluable deletion reads |
|---|---:|---:|---:|---:|---:|
| modify | 43,532 | 60 | 1,992,656 | 1,330,120 | 1,324,059 |
| inhibitor | 27,771 | 54 | 1,997,053 | 904,080 | 902,652 |

Original Excel filenames, SHA-256 hashes and run settings are in `data/source_metadata.json`. Excel formatting, unrelated sheets and original computer paths were not copied into these portable tables. These files are processed SIQ event data, not sequencing FASTQ files.

## Output folders

Each dataset gets `mh/`, `summary/` and `patterns/` under the chosen output directory. CSV files retain full numerical precision for downstream figure editing. Columns ending in `fraction` are proportions from 0 to 1; multiply by 100 for percent. Group summaries report arithmetic means, sample SDs and valid-sample counts. Blank/NA is distinct from numerical zero.

- `mh/events.csv`: row-level strict MH annotation with exclusion/review flags.
- `summary/`: sample/group event composition, pure-deletion length distributions and means, original-coordinate deletion coverage, strict MH distributions and length-specific MH fractions; modification comparisons also have descriptive length-standardized tables and weights.
- `patterns/pattern_definitions.csv`: reconstructed alleles, source/display coordinates, equivalent starts, legacy IDs and MH candidates.
- `patterns/all_*`: all-length deletion pattern counts, fractions, Top5 and coverage for the corresponding dataset.
- `patterns/length_*`: modification-only 6/12-bp versions, with exact-length denominators and separate strict-evaluable/QC counts.
- `figures/`: optional plots exported from the tables. Figure generation does not change analysis values.

Some large pattern tables contain only patterns observed somewhere within a condition/time group, with zero-filled sample counts within that group. They are not evidence that every unlisted allele is biologically impossible. The full run is suitable for re-analysis; `expected/` stores compact reference checks rather than another copy of every output.

## Changing inputs

For new SIQ workbooks, use `scripts/import_siq.py` to extract fields and validate reference-related settings. Check the exported sample manifest; alias conventions are study-specific. For a new reference, update `config/junction.json` and provide a matching catalog instead of applying historical allele IDs to unrelated sequences. Keep a copy of the configuration with each analysis run.
