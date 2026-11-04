# Validation of the prepared package

Validated on 2026-11-04 with Python 3.12.14. The complete `run_analysis.py --figures` workflow was executed from the bundled portable inputs for both datasets.

- **32 scientific unit tests passed**, including primer boundaries, exact Raw mismatches, internally located deletions, MH=0 versus NA, read weighting, equal-weight biological replicates, zero denominators, invalid coordinates and common-support standardization.
- **71,303 original rawData rows** were independently checked against the archived MH annotation. The portable annotation output matched all 30 archived event columns, allowing equivalent integer/float text formatting.
- **24,531 reference items matched** between the regenerated outputs and compact historical fixtures. Counts and identifiers were compared exactly; floating-point values used absolute tolerance `1e-10`.
- All-length Top5 reproduced **190 condition/time/rank rows** (100 modifications and 90 inhibitors).
- The 6/12-bp analysis reproduced **17 pattern definitions, 120 sample denominators, 340 group-pattern records and 194 positive Top5 selections**.
- Strict MH sample/group summaries, selected exact-length MH comparisons and descriptive length-standardization outputs matched the earlier analyses.
- Five static figures were generated as SVG, PDF and PNG. Representative layouts were visually inspected.

`expected/package_verification.json` records the reference comparison. `expected/independent_source_verification.json` records the original workbook-based check. These are preparation-time results; users should rerun the verification after changing code or inputs.

```bash
python3 scripts/run_analysis.py --outdir results
python3 scripts/verify_reference_results.py --results results
python3 -m unittest discover -s tests -v
```

Verification supports numerical reproducibility of the defined post-processing. It does not validate the underlying sequencing experiment or establish a biochemical repair mechanism. Inferential statistical analyses and final manuscript figure selection are outside this package's current scope.
