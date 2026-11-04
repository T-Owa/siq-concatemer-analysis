#!/usr/bin/env python3
"""Reproduce 20 boundary and interpretation checks for this donor assay.

Run at the repository root:
    python -m unittest discover -s tests -v

Only the Python standard library is needed for these checks. They create
synthetic deletion sequences from config.json and do not read or modify
the experimental Excel workbook. These examples validate the supplied
196 bp reference and J=136 configuration, not an arbitrary replacement.
"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
from analyze_junction_mh import annotate_deletion, prepare_model


class JunctionHomologyExamples(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config_path = REPOSITORY / "config" / "junction.json"
        config = json.loads(config_path.read_text(encoding="utf-8-sig"))
        cls.model = prepare_model(config)
        cls.reference = cls.model["ref"]
        actual = (
            len(cls.reference), cls.model["junction"],
            cls.model["left_primer_end"], cls.model["right_primer_start"],
        )
        if actual != (196, 136, 20, 176):
            raise ValueError(
                "These examples require the original 196 bp / J=136 "
                "configuration with primer boundaries 20 and 176."
            )

    def raw_for(self, start, end):
        return self.reference[:start] + self.reference[end:]

    def check_example(self, start, end, expected_h, expected_status,
                      expected_sequence=None, expected_bounds=None):
        result = annotate_deletion(
            self.raw_for(start, end), start, end, end - start, self.model,
        )
        self.assertEqual(result["analysis_status"], expected_status)
        if expected_h is None:
            self.assertFalse(result["MH_evaluable"])
            # An unevaluable event must remain missing, never become MH=0.
            self.assertIsNone(result.get("junction_MH_max"))
            return result

        self.assertTrue(result["MH_evaluable"])
        self.assertEqual(result["junction_MH_max"], expected_h)
        self.assertEqual(len(result["junction_MH_sequence"]), expected_h)
        a, b = result["candidate_start_min"], result["candidate_start_max"]
        d, j = end - start, self.model["junction"]
        self.assertEqual(b - a, expected_h)
        self.assertGreaterEqual(a, self.model["left_primer_end"])
        self.assertLessEqual(b + d, self.model["right_primer_start"])
        self.assertLessEqual(b, j)
        self.assertGreaterEqual(a + d, j)
        self.assertEqual(self.reference[a:b], self.reference[a + d:b + d])
        self.assertGreaterEqual(result["tail_loss_min"], 0)
        self.assertGreaterEqual(result["head_loss_min"], 0)
        if expected_sequence is not None:
            self.assertEqual(result["junction_MH_sequence"], expected_sequence)
        if expected_bounds is not None:
            self.assertEqual((a, b), expected_bounds)
        return result

    def test_01_head_internal_CCT_is_not_interdonor_MH3(self):
        """Both CCT copies are in the head: junction candidate MH=0."""
        result = self.check_example(
            139, 144, 0, "junction_compatible_with_internal_alternative",
            "", (136, 136),
        )
        self.assertEqual(result["equiv_start_count"], 4)
        self.assertTrue(result["internal_alternative_primer_intact"])

    def test_02_opposite_donor_GCC_is_MH3(self):
        """GCC occurs once on each donor side in the compatible overlap."""
        self.check_example(86, 174, 3, "junction_only", "GCC", (83, 86))

    def test_03_left_primer_boundary_is_inclusive(self):
        """Deleting immediately after the complete forward primer is allowed."""
        self.check_example(20, 136, 0, "junction_only", "", (20, 20))

    def test_04_right_primer_boundary_is_inclusive(self):
        """An exclusive deletion end at the reverse primer start is allowed."""
        self.check_example(136, 176, 0, "junction_only", "", (136, 136))

    def test_05_internal_only_deletion_is_unevaluable(self):
        self.check_example(145, 146, None, "no_single_junction_placement")

    def test_06_right_primer_overlap_without_rescue_is_unevaluable(self):
        """J-compatible placements exist, but none preserve both primers."""
        result = self.check_example(
            138, 180, None, "no_primer_preserving_junction_placement",
        )
        self.assertGreater(result["junction_candidate_count"], 0)
        self.assertEqual(result["valid_candidate_count"], 0)

    def test_07_one_base_inside_left_primer_is_unevaluable(self):
        self.check_example(
            19, 136, None, "no_primer_preserving_junction_placement",
        )

    def test_08_equivalent_placements_can_all_violate_left_primer(self):
        result = self.check_example(
            18, 137, None, "no_primer_preserving_junction_placement",
        )
        self.assertEqual(result["equiv_start_count"], 2)

    def test_09_primer_constraint_reduces_MH2_to_MH1(self):
        """Three J placements become two primer-preserving placements."""
        result = self.check_example(
            136, 177, 1, "junction_only", "G", (134, 135),
        )
        self.assertEqual(result["junction_candidate_count"], 3)
        self.assertEqual(result["valid_candidate_count"], 2)

    def test_10_tail_internal_GGT_is_not_interdonor_MH3(self):
        self.check_example(
            91, 136, 0, "junction_compatible_with_internal_alternative",
            "", (91, 91),
        )

    def test_11_CTCT_MH4_retains_left_right_assignment_ranges(self):
        result = self.check_example(
            111, 163, 4, "junction_only", "CTCT", (107, 111),
        )
        self.assertEqual((result["tail_loss_min"], result["tail_loss_max"]), (25, 29))
        self.assertEqual((result["head_loss_min"], result["head_loss_max"]), (23, 27))
        # Matching endpoints are paired; the two maxima are not additive.
        self.assertEqual(result["tail_loss_min"] + result["head_loss_max"], 52)
        self.assertEqual(result["tail_loss_max"] + result["head_loss_min"], 52)

    def test_12_single_terminal_base_deletion_is_MH0(self):
        self.check_example(135, 136, 0, "junction_only", "", (135, 135))

    def test_13_full_Raw_mismatch_is_not_silently_accepted(self):
        raw = self.raw_for(86, 174)
        changed_first_base = "C" if raw[0] != "C" else "A"
        result = annotate_deletion(
            changed_first_base + raw[1:], 86, 174, 88, self.model,
        )
        self.assertEqual(result["analysis_status"], "raw_reference_mismatch")
        self.assertFalse(result["MH_evaluable"])
        self.assertIsNone(result.get("junction_MH_max"))

    def test_14_inconsistent_deletion_length_is_rejected(self):
        result = annotate_deletion(self.raw_for(86, 174), 86, 174, 87, self.model)
        self.assertEqual(result["analysis_status"], "invalid_deletion_coordinates")
        self.assertFalse(result["MH_evaluable"])
        self.assertIsNone(result.get("junction_MH_max"))

    def test_15_fractional_coordinates_are_rejected(self):
        result = annotate_deletion(
            self.raw_for(86, 174), 86.5, 174, 87.5, self.model,
        )
        self.assertEqual(result["analysis_status"], "invalid_deletion_coordinates")
        self.assertFalse(result["MH_evaluable"])
        self.assertIsNone(result.get("junction_MH_max"))


    def test_16_six_bp_Cy5_pattern_has_MH2_TG(self):
        self.check_example(134, 140, 2, "junction_only", "TG", (132, 134))

    def test_17_twelve_bp_Cy5_pattern_has_MH3_CTG(self):
        self.check_example(128, 140, 3, "junction_only", "CTG", (125, 128))

    def test_18_D06_01_internal_tail_pattern_stays_NA(self):
        result = self.check_example(129, 135, None, "no_single_junction_placement")
        self.assertEqual((result["equiv_start_min"], result["equiv_start_max"]), (126, 129))

    def test_19_D06_06_internal_head_pattern_stays_NA(self):
        result = self.check_example(138, 144, None, "no_single_junction_placement")
        self.assertEqual(result["equiv_start_count"], 1)

    def test_20_truncated_terminal_A_is_not_silently_trimmed(self):
        # Mirrors the important QC distinction in the 6/12 bp Top5 population:
        # an SIQ-coordinate-reconstructed allele can be plotted, but a full Raw
        # with one extra terminal base missing cannot enter strict MH counts.
        raw = self.raw_for(134, 140)
        self.assertEqual(raw[-1], "A")
        result = annotate_deletion(raw[:-1], 134, 140, 6, self.model)
        self.assertEqual(result["analysis_status"], "raw_reference_mismatch")
        self.assertFalse(result["MH_evaluable"])
        self.assertIsNone(result.get("junction_MH_max"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
