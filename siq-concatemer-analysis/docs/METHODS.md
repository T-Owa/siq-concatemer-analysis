# Analysis definitions

## Reference and source settings

The ideal head-to-tail junction amplicon is 196 nt. The junction index is **J = 136** in 0-based coordinates: the last Tail base has index 135 and the first Head base index 136. Deletions use half-open intervals `[s,e)` with length `d = e - s`. The left primer occupies `[0,20)` and the reverse complement of the right primer occupies `[176,196)`. The complete reference and primers are in `config/junction.json`.

Both source workbooks report SIQ version 2.1, minimum support 2, maximum base error 0.08, TINS search distance 100, maximum CPUs 8, long-read analysis false, remerge false and filtering of DELINS caused by multiple events true. Per-sample `bases pastprimer` is 0. Source SIQ settings are recorded in `data/source_metadata.json`. Reported event counts have already undergone SIQ filtering; absent rows cannot be recovered by this code.

## Strict junction-compatible MH

Only positive-count pure `DELETION` rows are considered. Require integer, valid deletion coordinates, `d = e - s`, and exact equality of the entire observed Raw sequence Q with `R[:s] + R[e:]`. No trimming, mismatch tolerance or sequence correction is applied.

For fixed length d, enumerate all starts:

```text
S = {s : Q == R[:s] + R[s+d:]}
V = {s in S : 20 <= s <= 136 <= s+d <= 176}
```

The bounds permit one-sided deletions that touch the junction as well as deletions that remove sequence on both sides. If V is empty, the MH result is unevaluable (NA). Otherwise, with `a = min(V)` and `b = max(V)`, define:

```text
MH_max = b - a
MH_sequence = R[a:b] = R[a+d:b+d]
```

Both matching copies lie on opposite sides of the junction under these constraints. This restricts the annotation to junction-compatible homology, rather than treating every within-donor repeat ambiguity as junction MH. A repetitive allele may also have internal-deletion alternatives; those are retained as flags, not hidden.

The observed maximum in these data is 5 bp; the exact 5-bp category is not a `5+` bin. Threshold summaries explicitly use `>=1`, `>=2` or `>=3`. No upper MH cap of 5 is imposed by the algorithm.

The strict denominator is the count of evaluable pure deletions in each sample. MH=0 is included. Raw disagreements, invalid coordinates and lack of a primer-preserving junction-compatible placement are excluded from this denominator and reported separately. Sequence compatibility does not demonstrate biochemical MH use, MMEJ, TMEJ or a unique end-processing mechanism.

## Event, deletion length and position summaries

Event counts retain every Type reported by SIQ. The source workbooks contain WT, DELETION, INSERTION, SNV, DELINS, TINS, DELINS_SNV, TANDEMDUPLICATION and TANDEMDUPLICATION_COMPOUND. No type is silently merged into Other.

All-event counts are summed over `rawData`. They are not the input sequencing total or MergedCorrect. Pure deletion length distributions and weighted sample means use Type=DELETION. Strict-MH-evaluable deletion means are calculated separately and named separately.

Position coverage counts a deletion at base x if `s <= x < e`, using original SIQ coordinates. It is provided relative to pure deletion counts and all-event counts. It differs from an endpoint-frequency map. Repetitive sequences permit equivalent coordinate placements, so the coverage map is representation-dependent; it is not a map of uniquely observed biochemical cleavage positions.

## Deletion patterns and Top5

Pure deletion rows are grouped by identical coordinate-reconstructed reference alleles. Raw disagreements remain in the SIQ-based pattern counts and are flagged. Invalid-coordinate rows cannot be assigned a reconstructed allele; their counts remain in the applicable denominators and are reported as unclassified.

For the all-length analysis, rank each condition/time group's alleles by the arithmetic mean of the three sample-level within-pure-deletion fractions. The same selected patterns also receive all-event fractions. For exact-length analysis, select 6-bp or 12-bp deletions in the modification dataset and rank by fractions of that selected-length population. Only positive-frequency patterns can enter Top5. Groups with fewer than five observed patterns do not receive invented zero-frequency entries.

Ties use pattern ID order. Missing pattern counts in valid samples are zero; zero denominators yield NA. Top5 coverage is calculated by summing the same group-selected patterns within each sample, then computing the group mean and SD. Other retains the remaining fraction. Neither Top5 nor a selected display subset is renormalized to 100%.

Historical D06/D12 and DEL identifiers are linked through `config/pattern_catalog.csv`. IDs label alleles, not repair mechanisms. The source representative uses the maximum observed SIQ delStart for that reconstructed allele. Equivalent display coordinates are chosen by minimum distance of the interval to J, then minimum shift from this representative. Internal-only alleles remain internal; they are never moved to the junction if this would change the resulting sequence. The sequence-compatible MH annotation on reconstructed alleles is distinct from strict MH evaluation of full Raw.

In the archived 6/12-bp comparison, the two Raw disagreements retained for pattern counts are Acryd-12h-1 (6 bp, 663 reads) and sp18-12h-2 (12 bp, 107 reads). Both have a missing terminal A relative to the SIQ-coordinate reconstruction. These records are excluded from strict MH evaluation. D06-01 and D06-06 have exact Raw matches but no junction-compatible deletion placement; their MH is NA.

## Replicates and descriptive standardization

For each condition/time, the supplied data contain three independent biological experiments. Sample labels with matching replicate numbers do not imply pairing between conditions. Counts are weights within a sample; individual sequencing reads are not independent biological replicates. Group SD uses `ddof=1`. The number of valid denominators is reported for each metric.

For a descriptive Cy5 comparison at a given time, use exact deletion lengths having at least 20 strict-MH-evaluable reads in every sample of both conditions (six samples here). Normalize each sample's length counts over these common lengths. The mean of the six normalized length distributions defines the shared weights. Apply these weights to each sample's length-specific MH fraction, then compare equal-weight condition means. Coverage of the retained common lengths is reported. Optional thresholds of 1 or 100 reads can be requested explicitly for sensitivity analysis.

Matching total deletion length does not match left/right allocation or endpoint sequence preference. Residual differences are not proof of a repair-pathway switch. No p-values or significance claims are generated by this package.
