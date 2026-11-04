# Historical reference results

These compact CSV fixtures were extracted from the earlier analysis deliverables, independently of the new publication-package aggregation code. Their source archive/workbook hashes are recorded in `fixture_provenance.json`.

- `historical_mh_summary_*`: original strict-MH per-sample summaries (60 modification and 54 inhibitor samples).
- `historical_sample_metrics.csv` and `historical_group_metrics.csv`: original strict-MH length/MH summaries and pure-DELETION frequencies.
- `historical_exact_length_mh_samples.csv`: original sample-level strict-MH results at 6, 12, 17, 19, 64 and 88 bp for both datasets.
- `historical_standardized.csv`: original Cy5 comparisons at a minimum of 20 reads per length in every sample. The original percent scale is retained; the verification script explicitly converts it to fractions.
- `historical_all_top5.csv`: 190 Top5 rows from the modification/inhibitor workbook, retaining historical DEL IDs and source-coordinate labels.
- `historical_length_*`: original 6/12-bp modification-only pattern definitions, sample denominators, group frequencies and coverage. Internal-only and Raw-disagreement reads remain in those denominators.
- `independent_source_verification.json`: separate direct re-enumeration against both original SIQ workbooks. All 71,303 rawData rows and the archived 6/12-bp statistics were checked. No analysis implementation was imported for that verification.

Run from the repository root after producing results:

```bash
python3 scripts/verify_reference_results.py --results results --report results/verification.json
```

The verification compares scientific values, with exact integer counts and a floating-point tolerance of 1e-10. It ignores timestamps, file paths and cosmetic labels. Equivalent display coordinates may differ from the old workbook's source coordinates; the reconstructed allele, historical ID and quantitative results are checked. This verification applies to the supplied study data, not arbitrary replacement data.
