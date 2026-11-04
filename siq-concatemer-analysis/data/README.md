# Analysis inputs

These are minimal exports of the **SIQ `rawData` worksheet**, not raw FASTQ
files. They retain all reported rows, including zero-count reference rows.
No additional event-support threshold was applied during export. Events
omitted by SIQ cannot be recovered from these tables.

| File | Dataset | Samples | SIQ rows | Sum of `countEvents` |
|---|---|---:|---:|---:|
| `modify_events.csv.gz` | Donor modifications | 60 | 43,532 | 1,992,656 |
| `inhibitor_events.csv.gz` | Repair inhibitors | 54 | 27,771 | 1,997,053 |

The CSV files are UTF-8 with a BOM and deterministic gzip compression.
Python scripts read them directly without manual decompression.

## Event table fields

| Field | Meaning |
|---|---|
| `source_excel_row` | Original 1-based worksheet row number; header is row 1 |
| `Alias` | Original sample identifier |
| `Type` | Original SIQ event classification |
| `countEvents` | Number of SIQ-reported reads supporting the row |
| `Raw` | Complete SIQ sequence, normalized to uppercase |
| `delStart`, `delEnd` | Original SIQ deletion coordinates, interpreted as a 0-based half-open interval |
| `delSize` | SIQ deletion size in bp |
| `homology`, `homologyLength` | Original SIQ annotations; these are distinct from the custom junction-MH annotation |

`samples.csv` maps each Alias to dataset, condition, display label, collection
time and biological replicate. Replicate numbers identify independent
experiments and do **not** imply pairing between conditions. Original aliases
are retained; labels such as `Bractppin` are normalized only in
`condition_label`.

`source_metadata.json` records the original workbook filenames and SHA-256
hashes, export hashes, checked primer/flank settings and SIQ run settings.
Unrelated workbook columns and original filesystem paths are omitted.

The core MH annotation requires the complete observed `Raw` to match the
SIQ-coordinate reconstruction. The pattern analysis also retains
coordinate-reconstructable events whose `Raw` does not match, with explicit
QC fields. See the main methods documentation for the distinct denominators.
