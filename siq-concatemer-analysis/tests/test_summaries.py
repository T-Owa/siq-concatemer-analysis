import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from summarize_events import fraction, group_rows, read_samples, standardize, summarize
from analyze_junction_mh import prepare_model


class SummaryTests(unittest.TestCase):
    def setUp(self):
        self.model = prepare_model(json.loads((ROOT / "config/junction.json").read_text()))
        self.samples = [dict(dataset="modify", Alias=f"X-{i}", condition="X", condition_label="X", time="12h", replicate=i)
                        for i in [1, 2, 3]]

    def analyze(self, rows, samples=None):
        fields = ["Alias", "Type", "countEvents", "delStart", "delEnd", "delSize", "MH_evaluable",
                  "junction_MH_max", "analysis_status"]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "events.csv"
            with path.open("w", newline="") as h:
                writer = csv.DictWriter(h, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            return summarize(path, samples or self.samples, self.model)

    def test_populations_and_half_open_positions(self):
        rows = [dict(Alias="X-1", Type="WT", countEvents=10, MH_evaluable="False"),
                dict(Alias="X-1", Type="DELETION", countEvents=2, delStart=134, delEnd=140, delSize=6,
                     MH_evaluable="True", junction_MH_max=2, analysis_status="junction_only"),
                dict(Alias="X-1", Type="DELETION", countEvents=8, delStart=129, delEnd=135, delSize=6,
                     MH_evaluable="False", analysis_status="no_single_junction_placement"),
                dict(Alias="X-2", Type="WT", countEvents=100, MH_evaluable="False")]
        tables, state, count = self.analyze(rows)
        sample = tables["sample_metrics"][0]
        self.assertEqual(count, 4)
        self.assertEqual(sample["pure_deletion_reads"], 10)
        self.assertEqual(sample["evaluable_deletion_reads"], 2)
        self.assertEqual(sample["unevaluable_deletion_reads"], 8)
        self.assertEqual(sample["deletion_fraction_of_reported"], .5)
        self.assertEqual(sample["mh_ge2_fraction_of_evaluable"], 1)
        self.assertEqual(state["X-1"]["positions"][134], 10)
        self.assertEqual(state["X-1"]["positions"][135], 2)
        self.assertEqual(state["X-1"]["positions"][140], 0)
        group = next(r for r in tables["event_types_by_group"] if r["Type"] == "DELETION")
        self.assertEqual(group["event_fraction_mean"], .25)  # sample 3 has no reported reads: NA
        self.assertEqual(group["event_fraction_n_valid"], 2)
        self.assertEqual(group["n_samples"], 3)
        absent_length = next(r for r in tables["exact_length_mh_by_sample"]
                             if r["Alias"] == "X-1" and r["deletion_length_bp"] == 12)
        self.assertIsNone(absent_length["mh_ge2_fraction"])

    def test_invalid_coordinates_remain_in_deletion_denominator(self):
        tables, state, _ = self.analyze([dict(Alias="X-1", Type="DELETION", countEvents=3,
                                    delStart=0, delEnd=999, delSize=999, MH_evaluable="False",
                                    analysis_status="invalid_deletion_coordinates")])
        record = tables["sample_metrics"][0]
        self.assertEqual(record["pure_deletion_reads"], 3)
        self.assertEqual(record["invalid_length_deletion_reads"], 3)
        self.assertIsNone(record["mean_deletion_bp_all_deletions"])
        self.assertEqual(sum(state["X-1"]["positions"].values()), 0)

    def test_unknown_alias_and_negative_counts_rejected(self):
        for row in [dict(Alias="missing", Type="WT", countEvents=1, MH_evaluable="False"),
                    dict(Alias="X-1", Type="WT", countEvents=-1, MH_evaluable="False")]:
            with self.assertRaises(ValueError):
                self.analyze([row])

    def test_read_weighted_mean_and_equal_experiment_average(self):
        rows = [dict(Alias="X-1", Type="DELETION", countEvents=9, delStart=135, delEnd=136, delSize=1,
                     MH_evaluable="True", junction_MH_max=0, analysis_status="junction_only"),
                dict(Alias="X-1", Type="DELETION", countEvents=1, delStart=130, delEnd=140, delSize=10,
                     MH_evaluable="True", junction_MH_max=0, analysis_status="junction_only"),
                dict(Alias="X-2", Type="DELETION", countEvents=1000, delStart=134, delEnd=136, delSize=2,
                     MH_evaluable="True", junction_MH_max=0, analysis_status="junction_only")]
        tables, _, _ = self.analyze(rows)
        self.assertAlmostEqual(tables["sample_metrics"][0]["mean_deletion_bp_all_deletions"], 1.9)
        self.assertAlmostEqual(tables["group_metrics"][0]["mean_deletion_bp_all_deletions_mean"], 1.95)

    def test_standardization_common_support_and_weights(self):
        samples, state = [], {}
        for condition in ["Cy5", "X"]:
            for replicate in [1, 2, 3]:
                alias = f"{condition}-{replicate}"
                samples.append(dict(dataset="modify", Alias=alias, condition=condition, condition_label=condition,
                                    time="12h", replicate=replicate))
                exact = {6: {0: 30, 2: 10}, 12: {0: 10, 2: 10}} if condition == "Cy5" else {6: {0: 10, 2: 10}, 12: {0: 40}}
                exact[17] = {2: 1}
                state[alias] = {"exact": exact, "strict": 61}
        result, individual, weights = standardize(samples, state, [20])
        row = next(r for r in result if r["mh_threshold"] == 2)
        self.assertEqual(row["common_lengths"], "6;12")
        self.assertAlmostEqual(row["focus_adjusted_fraction_mean"], .375)
        self.assertAlmostEqual(row["comparison_adjusted_fraction_mean"], .25)
        self.assertAlmostEqual(row["difference_fraction"], .125)
        self.assertAlmostEqual(row["focus_coverage_fraction_mean"], 60 / 61)
        self.assertEqual([w["common_weight_fraction"] for w in weights], [.5, .5])

    def test_no_common_lengths_is_reported_as_na(self):
        samples = [dict(dataset="modify", Alias=c, condition=c, condition_label=c, time="12h", replicate=1)
                   for c in ["Cy5", "X"]]
        state = {"Cy5": {"exact": {6: {2: 10}}, "strict": 10}, "X": {"exact": {12: {3: 10}}, "strict": 10}}
        results, _, weights = standardize(samples, state, [20])
        self.assertEqual(weights, [])
        self.assertEqual(results[0]["status"], "no_common_lengths")
        self.assertIsNone(results[0]["focus_adjusted_fraction_mean"])
        self.assertIsNone(fraction(0, 0))


if __name__ == "__main__":
    unittest.main()
