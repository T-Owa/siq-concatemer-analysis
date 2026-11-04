#!/usr/bin/env python3
"""Read-only SIQ workbook analysis of maximum junction-compatible homology.

This is a custom sequence-based annotation, not a SIQ option and not proof
of MH use or of a particular repair pathway. See README_JA.txt.
Python >=3.10; only external dependency: openpyxl.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import shutil
import sys
import tempfile
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VERSION = "1.0.0"
EVENT_FIELDS = [
    "source_excel_row", "Alias", "Type", "countEvents", "Raw",
    "delStart", "delEnd", "delSize", "siq_homology", "siq_homologyLength",
    "analysis_status", "MH_evaluable", "junction_MH_max", "junction_MH_sequence",
    "raw_exact_single_deletion", "equiv_start_min", "equiv_start_max",
    "equiv_start_count", "junction_candidate_count", "valid_candidate_count",
    "candidate_start_min", "candidate_start_max",
    "internal_alternative_any", "internal_alternative_primer_intact",
    "primer_intact_possible", "tail_loss_min", "tail_loss_max",
    "head_loss_min", "head_loss_max", "review_flag",
]


def dna(value):
    return str(value or "").strip().upper()


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def integer(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name}: expected an integer, got {value!r}") from None
    if not math.isfinite(number) or number != int(number):
        raise ValueError(f"{name}: expected an integer, got {value!r}")
    return int(number)


def unique_position(ref, seq, name):
    positions = [i for i in range(len(ref) - len(seq) + 1)
                 if ref[i:i + len(seq)] == seq]
    if not seq or len(positions) != 1:
        raise ValueError(f"{name} must occur exactly once in the reference: {positions}")
    return positions[0]


def prepare_model(config):
    ref = dna(config["reference_sequence"])
    left = dna(config["left_primer"])
    right = dna(config["right_primer"])
    if any(set(s) - set("ACGT") for s in (ref, left, right)):
        raise ValueError("Reference and primers must contain only A/C/G/T.")
    junction = integer(config["junction_index"], "junction_index")
    left_start = unique_position(ref, left, "left_primer")
    right_start = unique_position(ref, reverse_complement(right), "right_primer RC")
    left_end = left_start + len(left)
    if not (left_end <= junction <= right_start):
        raise ValueError("The junction must lie between the two primer-binding intervals.")
    return {"ref": ref, "junction": junction, "left_primer": left,
            "right_primer": right, "left_primer_start": left_start,
            "left_primer_end": left_end, "right_primer_start": right_start,
            "right_primer_end": right_start + len(right)}


def annotate_deletion(raw, start, end, size, model):
    """Enumerate all exact single-deletion placements; retain NA separately.

    R is the ideal joined reference and Q the full SIQ Raw string.
    For fixed d, S={s: Q == R[:s] + R[s+d:]}. Valid placements also satisfy
    left_primer_end <= s <= J <= s+d <= right_primer_start.
    If a=min(valid), b=max(valid), h=b-a and R[a:b]==R[a+d:b+d].
    The two copies lie on opposite sides of J. h is a maximum compatible
    sequence homology, not a measurement of actual biochemical MH use.
    """
    result = {"MH_evaluable": False, "analysis_status": "invalid_deletion_coordinates",
              "review_flag": "invalid_deletion_coordinates"}
    ref, junction = model["ref"], model["junction"]
    try:
        s = integer(start, "delStart")
        e = integer(end, "delEnd")
        d = integer(size, "delSize")
    except ValueError:
        return result
    if not (0 <= s < e <= len(ref) and d == e - s):
        return result
    raw = dna(raw)
    canonical_exact = ref[:s] + ref[e:] == raw
    result["raw_exact_single_deletion"] = canonical_exact
    if not canonical_exact:
        # Do not silently trim reads, tolerate mismatches, or infer MH from
        # SIQ coordinates alone when the complete observed Raw disagrees.
        result.update(analysis_status="raw_reference_mismatch",
                      review_flag="full_Raw_not_reconstructed_by_SIQ_deletion")
        return result
    starts = [i for i in range(len(ref) - d + 1)
              if ref[:i] + ref[i + d:] == raw]
    if starts != list(range(starts[0], starts[-1] + 1)):
        raise AssertionError("Non-contiguous placements violate the exact single-deletion model.")
    junction_starts = [i for i in starts if i <= junction <= i + d]
    primer_starts = [i for i in starts
                     if model["left_primer_end"] <= i
                     and i + d <= model["right_primer_start"]]
    valid = [i for i in primer_starts if i <= junction <= i + d]
    internal_any = any(not (i <= junction <= i + d) for i in starts)
    internal_intact = any(not (i <= junction <= i + d) for i in primer_starts)
    result.update(equiv_start_min=starts[0], equiv_start_max=starts[-1],
                  equiv_start_count=len(starts), junction_candidate_count=len(junction_starts),
                  valid_candidate_count=len(valid), internal_alternative_any=internal_any,
                  internal_alternative_primer_intact=internal_intact,
                  primer_intact_possible=bool(primer_starts))
    if not valid:
        reason = ("no_single_junction_placement" if not junction_starts
                  else "no_primer_preserving_junction_placement")
        result.update(analysis_status=reason, review_flag=reason)
        return result
    a, b = valid[0], valid[-1]
    h = b - a
    sequence = ref[a:b]
    assert sequence == ref[a + d:b + d]
    assert b <= junction <= a + d and 0 <= h <= d
    result.update(
        MH_evaluable=True,
        analysis_status=("junction_compatible_with_internal_alternative"
                         if internal_any else "junction_only"),
        junction_MH_max=h, junction_MH_sequence=sequence,
        candidate_start_min=a, candidate_start_max=b,
        tail_loss_min=junction - b, tail_loss_max=junction - a,
        head_loss_min=a + d - junction, head_loss_max=b + d - junction,
        review_flag=("also_explainable_by_internal_deletion" if internal_any else ""),
    )
    return result


def worksheet_records(ws):
    # Some SIQ workbooks have incomplete XML dimension metadata.
    ws.reset_dimensions()
    iterator = ws.iter_rows(values_only=True)
    try:
        headers = [str(v).strip() if v is not None else "" for v in next(iterator)]
    except StopIteration:
        raise ValueError(f"Empty sheet: {ws.title}") from None
    for excel_row, values in enumerate(iterator, start=2):
        if all(v is None for v in values):
            continue
        yield excel_row, dict(zip(headers, values))


def validate_runinfo(workbook, config, model):
    if "RunInfo" not in workbook.sheetnames:
        raise ValueError("RunInfo sheet is required to verify primers and junction settings.")
    result = []
    for rownum, row in worksheet_records(workbook["RunInfo"]):
        alias = str(row.get("alias") or "").strip()
        if not alias:
            raise ValueError(f"RunInfo row {rownum}: missing alias")
        for column, value in (("left primer", model["left_primer"]),
                              ("right primer", model["right_primer"])):
            if dna(row.get(column)) != value:
                raise ValueError(f"{alias}: {column} differs from config.json")
        left_flank, right_flank = dna(row.get("left flank")), dna(row.get("right flank"))
        j, ref = model["junction"], model["ref"]
        if not left_flank or not ref[:j].endswith(left_flank):
            raise ValueError(f"{alias}: left flank does not end at the configured junction")
        if not right_flank or not ref[j:].startswith(right_flank):
            raise ValueError(f"{alias}: right flank does not start at the configured junction")
        past = integer(row.get("#bases pastprimer"), "#bases pastprimer")
        if past != integer(config["expected_bases_pastprimer"], "expected_bases_pastprimer"):
            raise ValueError(f"{alias}: bases pastprimer={past}, differs from config.json")
        result.append({"Alias": alias, "bases_pastprimer": past,
                       "left_flank": left_flank, "right_flank": right_flank})
    if not result or len({r["Alias"] for r in result}) != len(result):
        raise ValueError("RunInfo aliases must be nonempty and unique.")
    return result


def write_csv(path, fields, records):
    # BOM allows Excel on Windows to recognize UTF-8 automatically.
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def run(input_path, config_path, outdir):
    import openpyxl

    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    model = prepare_model(config)
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise ValueError("Output folder is not empty. Choose a new --outdir; existing files are not overwritten.")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", message="Workbook contains no default style")
    workbook = openpyxl.load_workbook(input_path, read_only=True, data_only=True)
    try:
        runinfo = validate_runinfo(workbook, config, model)
        allowed_aliases = {r["Alias"] for r in runinfo}
        if "rawData" not in workbook.sheetnames:
            raise ValueError("No rawData sheet in the input workbook.")
        stats = {r["Alias"]: Counter() for r in runinfo}
        distributions = {r["Alias"]: Counter() for r in runinfo}
        rows_seen = 0
        required = {"countEvents", "Alias", "Type", "Raw", "delStart", "delEnd", "delSize"}
        with tempfile.TemporaryDirectory(prefix="mh_work_", dir=outdir.parent) as tmp:
            staging = Path(tmp)
            with (staging / "events.csv").open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=EVENT_FIELDS, extrasaction="ignore")
                writer.writeheader()
                for rownum, row in worksheet_records(workbook["rawData"]):
                    if not required <= row.keys():
                        raise ValueError(f"Missing required rawData columns: {sorted(required - row.keys())}")
                    alias = str(row["Alias"] or "").strip()
                    if alias not in allowed_aliases:
                        raise ValueError(f"rawData row {rownum}: Alias {alias!r} not in RunInfo")
                    count = integer(row["countEvents"], f"row {rownum} countEvents")
                    if count < 0:
                        raise ValueError(f"row {rownum}: negative read count")
                    kind, raw = str(row["Type"] or "").strip(), dna(row["Raw"])
                    if kind == "WT" and raw != model["ref"]:
                        raise ValueError(f"row {rownum}: WT Raw differs from the configured reference")
                    event = {k: row.get(k) for k in ("delStart", "delEnd", "delSize")}
                    event.update(source_excel_row=rownum, Alias=alias, Type=kind,
                                 countEvents=count, Raw=raw, siq_homology=row.get("homology"),
                                 siq_homologyLength=row.get("homologyLength"), MH_evaluable=False)
                    stats[alias]["total_reported_reads"] += count
                    if count == 0:
                        event["analysis_status"] = "zero_count_reference"
                    elif kind != "DELETION":
                        event["analysis_status"] = "not_target_type"
                    else:
                        stats[alias]["pure_deletion_reads"] += count
                        event.update(annotate_deletion(raw, row["delStart"], row["delEnd"],
                                                       row["delSize"], model))
                        if event["MH_evaluable"]:
                            h = event["junction_MH_max"]
                            stats[alias]["evaluable_deletion_reads"] += count
                            stats[alias]["MH_weighted_sum"] += h * count
                            distributions[alias][h] += count
                            if event["internal_alternative_any"]:
                                stats[alias]["internal_alternative_reads"] += count
                            if event["internal_alternative_primer_intact"]:
                                stats[alias]["internal_alternative_primer_intact_reads"] += count
                        else:
                            stats[alias]["excluded_" + event["analysis_status"] + "_reads"] += count
                    writer.writerow(event)
                    rows_seen += 1
                    if rows_seen % 10000 == 0:
                        print(f"Processed {rows_seen:,} rawData rows...", flush=True)
            if not rows_seen:
                raise ValueError("No rawData rows.")
            summary, distribution_rows = [], []
            maximum_h = max((max(d, default=0) for d in distributions.values()), default=0)
            excluded_keys = ["raw_reference_mismatch", "invalid_deletion_coordinates",
                             "no_single_junction_placement", "no_primer_preserving_junction_placement"]
            for alias, counters in stats.items():
                n, total = counters["evaluable_deletion_reads"], counters["pure_deletion_reads"]
                mh = distributions[alias]
                rec = {"Alias": alias,
                       "total_reported_reads": counters["total_reported_reads"],
                       "pure_deletion_reads": total, "evaluable_deletion_reads": n,
                       "unevaluable_deletion_reads": total - n,
                       "evaluable_fraction_of_deletions": n / total if total else None,
                       "MH0_reads": mh[0], "MH1_reads": mh[1],
                       "MH2to3_reads": mh[2] + mh[3],
                       "MHge4_reads": sum(v for k, v in mh.items() if k >= 4),
                       "MHge2_reads": sum(v for k, v in mh.items() if k >= 2),
                       "MHge2_fraction_of_evaluable": sum(v for k, v in mh.items() if k >= 2) / n if n else None,
                       "mean_junction_MH_bp": counters["MH_weighted_sum"] / n if n else None,
                       "internal_alternative_reads": counters["internal_alternative_reads"],
                       "internal_alternative_fraction_of_evaluable": counters["internal_alternative_reads"] / n if n else None,
                       "internal_alternative_primer_intact_reads": counters["internal_alternative_primer_intact_reads"]}
                for key in excluded_keys:
                    rec["excluded_" + key + "_reads"] = counters["excluded_" + key + "_reads"]
                assert sum(mh.values()) == n
                assert sum(rec["excluded_" + k + "_reads"] for k in excluded_keys) == total - n
                summary.append(rec)
                for h in range(maximum_h + 1):
                    distribution_rows.append({"Alias": alias, "junction_MH_bp": h,
                                              "read_count": mh[h], "evaluable_deletion_reads": n,
                                              "fraction_of_evaluable": mh[h] / n if n else None})
            write_csv(staging / "sample_summary.csv", list(summary[0]), summary)
            write_csv(staging / "mh_distribution.csv", list(distribution_rows[0]), distribution_rows)
            settings = []
            if "SIQRun" in workbook.sheetnames:
                settings = [row for _, row in worksheet_records(workbook["SIQRun"])]
            metadata = {
                "script_version": VERSION, "utc_time": datetime.now(timezone.utc).isoformat(),
                "python_version": platform.python_version(), "openpyxl_version": openpyxl.__version__,
                "input_filename": input_path.name, "input_sha256": digest(input_path),
                "script_sha256": digest(Path(__file__)), "config_sha256": digest(config_path),
                "config": config, "model": model, "SIQRun": settings,
                "runinfo_checked": runinfo, "rawData_rows": rows_seen,
                "n_samples": len(stats),
                "total_reported_reads": sum(x["total_reported_reads"] for x in summary),
                "pure_deletion_reads": sum(x["pure_deletion_reads"] for x in summary),
                "evaluable_deletion_reads": sum(x["evaluable_deletion_reads"] for x in summary),
                "unevaluable_deletion_reads": sum(x["unevaluable_deletion_reads"] for x in summary),
                "coordinate_convention": f"0-based half-open [s,e); J={model['junction']} means bases {model['junction']}|{model['junction'] + 1} in 1-based coordinates",
                "MH_definition": "Maximum exact homology compatible with the observed full Raw, a single fixed-length deletion across the configured donor junction, and preservation of both primer-binding intervals.",
                "denominator": "Only positive-count DELETION reads with MH_evaluable=True; MH=0 included; missing/unevaluable values excluded and reported separately.",
                "interpretation": "Sequence-compatible maximum only, not proof of actual MH use, MMEJ, or TMEJ; internal-deletion alternatives retained as flags.",
                "no_extra_support_filter": "Uses countEvents already reported in rawData; omitted SIQ events cannot be recovered.",
                "sources": ["https://github.com/RobinVanSchendel/SIQ", "https://doi.org/10.1093/nargab/lqac063"],
            }
            (staging / "run_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
            if outdir.exists() and any(outdir.iterdir()):
                raise ValueError("Output folder became nonempty during analysis; nothing was overwritten.")
            outdir.mkdir(exist_ok=True)
            for path in staging.iterdir():
                shutil.move(str(path), str(outdir / path.name))
            print(f"Done: {len(stats)} samples; {metadata['evaluable_deletion_reads']:,} evaluable deletion reads.")
            print(f"Unevaluable deletion reads: {metadata['unevaluable_deletion_reads']:,}")
            print(f"Results: {outdir}")
    finally:
        workbook.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Original SIQ .xlsx workbook (not modified)")
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config.json"))
    parser.add_argument("--outdir", type=Path, default=Path("mh_results"), help="New or empty output directory")
    args = parser.parse_args()
    if sys.version_info < (3, 10):
        parser.error("Python 3.10 or newer is required.")
    try:
        run(args.input.resolve(), args.config.resolve(), args.outdir.resolve())
    except (ValueError, OSError, ImportError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
