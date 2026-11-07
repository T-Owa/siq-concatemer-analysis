"""Numerical and denominator tests for reconstructed deletion-pattern summaries."""
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from aggregate_patterns import aggregate, load_csv, mean_sd

CONFIG = json.loads((ROOT / "config/junction.json").read_text())
CATALOG = load_csv(ROOT / "config/pattern_catalog.csv")
REF = CONFIG["reference_sequence"]


def samples(n=3):
    return [dict(dataset="modify", Alias=f"test-24h-{i}", condition="test",
                 condition_label="Test", time="24h", replicate=str(i))
            for i in range(1, n + 1)]


def event(rep, start, length, count, raw=None):
    end = start + length
    return dict(Alias=f"test-24h-{rep}", Type="DELETION", countEvents=str(count),
                delStart=str(start), delEnd=str(end), delSize=str(length),
                Raw=REF[:start] + REF[end:] if raw is None else raw)


def wt(rep, count=1):
    return dict(Alias=f"test-24h-{rep}", Type="WT", countEvents=str(count),
                delStart="", delEnd="", delSize="", Raw=REF)


def run(events):
    return aggregate(events, samples(), CATALOG, CONFIG, "modify")


class PatternTests(unittest.TestCase):
    def test_equal_weight_samples_not_pooled_reads(self):
        rows = [event(1, 134, 6, 9), event(1, 136, 6, 1),
                event(2, 134, 6, 10), event(2, 136, 6, 90),
                event(3, 136, 6, 1), wt(1, 90)]
        result = run(rows)
        a = next(r for r in result["length_group_patterns"] if r["pattern_id"] == "D06-03")
        self.assertAlmostEqual(a["mean_fraction"], 1/3)
        self.assertNotAlmostEqual(a["mean_fraction"], 19/111)
        self.assertEqual(a["replicate_counts"], [9, 10, 0])
        self.assertAlmostEqual(a["mean_fraction_all_events"], (.09 + .1 + 0)/3)
        self.assertEqual(a["n_valid"], 3)
        self.assertAlmostEqual(a["sd_fraction"], mean_sd([.9, .1, 0])[1])

    def test_raw_mismatch_kept_in_pattern_denominator_not_strict_mh(self):
        raw = (REF[:134] + REF[140:])[:-1]
        result = run([event(1, 134, 6, 4, raw), event(2, 134, 6, 5), wt(3)])
        definition = next(r for r in result["pattern_definitions"] if r["pattern_id"] == "D06-03")
        self.assertEqual(definition["mh_length"], 2)
        self.assertEqual(definition["raw_mismatch_count"], 4)
        total = next(r for r in result["length_sample_totals"] if r["alias"] == "test-24h-1" and r["length"] == 6)
        self.assertEqual(total["denominator"], 4)
        self.assertEqual(total["strict_evaluable_denominator"], 0)
        self.assertEqual(total["raw_mismatch_count"], 4)

    def test_internal_allele_mh_is_missing_not_zero(self):
        result = run([event(1, 129, 6, 2), wt(2), wt(3)])
        d = result["pattern_definitions"][0]
        self.assertEqual(d["pattern_id"], "D06-01")
        self.assertIsNone(d["mh_length"])
        self.assertFalse(d["junction_compatible"])
        self.assertEqual(d["display_start_relative"], -7)
        self.assertEqual(d["display_end_relative"], -1)

    def test_zero_denominator_excluded_from_mean_and_rank(self):
        result = run([event(1, 134, 6, 1), event(2, 136, 6, 1), wt(3)])
        a = next(r for r in result["length_group_patterns"] if r["pattern_id"] == "D06-03")
        self.assertEqual(a["replicate_fractions"], [1., 0., None])
        self.assertEqual(a["n_valid"], 2)
        self.assertAlmostEqual(a["mean_fraction"], .5)
        coverage = next(r for r in result["length_coverage"] if r["length"] == 6)
        self.assertEqual(coverage["n_top_patterns"], 2)
        self.assertEqual(coverage["top5_fractions"], [1., 1., None])
        twelve = next(r for r in result["length_coverage"] if r["length"] == 12)
        self.assertIsNone(twelve["top5_mean_fraction"])
        self.assertEqual(twelve["n_valid"], 0)

    def test_invalid_coordinates_retained_in_denominator_and_other(self):
        invalid = event(1, 134, 6, 3)
        invalid["delStart"] = "invalid"
        result = run([event(1, 134, 6, 1), invalid, wt(2), wt(3)])
        total = next(r for r in result["length_sample_totals"] if r["alias"] == "test-24h-1" and r["length"] == 6)
        self.assertEqual(total["denominator"], 4)
        self.assertEqual(total["unclassified_count"], 3)
        coverage = next(r for r in result["length_coverage"] if r["length"] == 6)
        self.assertEqual(coverage["top5_mean_fraction"], .25)
        self.assertEqual(coverage["other_mean_fraction"], .75)

    def test_duplicate_sample_rejected(self):
        metadata = samples()
        metadata.append(metadata[0])
        with self.assertRaisesRegex(ValueError, "Duplicate sample"):
            aggregate([wt(1)], metadata, CATALOG, CONFIG, "modify")


if __name__ == "__main__":
    unittest.main()
