#!/usr/bin/env python3
"""Descriptive sample-level summaries of the strict-MH annotated SIQ events.

Rates are fractions (0..1); experiments, not reads, are the group-level units.
The all-DELETION and strict-MH populations are deliberately kept separate.
No hypothesis tests, pathway assignment, or amplicon ascertainment correction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

from analyze_junction_mh import integer, prepare_model

VERSION = "1.0.0"
SAMPLE_KEYS = ["dataset", "Alias", "condition", "condition_label", "time", "replicate"]
GROUP_KEYS = ["dataset", "condition", "condition_label", "time"]
THRESHOLDS = (1, 2, 3)


def fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def summary(values):
    values = [v for v in values if v is not None and math.isfinite(v)]
    return (statistics.mean(values) if values else None,
            statistics.stdev(values) if len(values) > 1 else None, len(values))


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_samples(path, dataset):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not set(SAMPLE_KEYS) <= set(reader.fieldnames or []):
            raise ValueError("Sample table requires: " + ", ".join(SAMPLE_KEYS))
        samples = []
        for row in reader:
            if row["dataset"] != dataset:
                continue
            sample = {k: str(row[k]).strip() for k in SAMPLE_KEYS}
            sample["replicate"] = integer(sample["replicate"], "replicate")
            if any(not sample[k] for k in SAMPLE_KEYS if k != "replicate"):
                raise ValueError("Sample metadata contains a missing field")
            samples.append(sample)
    if not samples:
        raise ValueError(f"No sample metadata for {dataset!r}")
    if len({s["Alias"] for s in samples}) != len(samples):
        raise ValueError("Duplicate sample Alias")
    if len({(s['condition'], s['time'], s['replicate']) for s in samples}) != len(samples):
        raise ValueError("Duplicate condition/time/replicate in sample table")
    for condition in {s["condition"] for s in samples}:
        if len({s["condition_label"] for s in samples if s["condition"] == condition}) != 1:
            raise ValueError(f"Inconsistent condition_label for {condition!r}")
    return samples


def group_rows(rows, extra_keys, metrics):
    """Equal-weight arithmetic mean and sample SD; missing denominators stay NA."""
    grouped = defaultdict(list)
    keys = GROUP_KEYS + list(extra_keys)
    for row in rows:
        grouped[tuple(row[k] for k in keys)].append(row)
    result = []
    for key, members in grouped.items():
        item = dict(zip(keys, key))
        item["n_samples"] = len(members)
        for metric in metrics:
            mean, sd, n = summary([r[metric] for r in members])
            item.update({metric + "_mean": mean, metric + "_sd": sd,
                         metric + "_n_valid": n})
        result.append(item)
    return result


def summarize(events_path, samples, model):
    """Aggregate counts; pure-DELETION includes internal and Raw-mismatch events."""
    state = {s["Alias"]: {
        "total": 0, "deletion": 0, "length_valid": 0, "length_sum": 0,
        "coord_invalid": 0, "strict": 0, "strict_length_sum": 0, "mh_sum": 0,
        "types": Counter(), "lengths": Counter(), "positions": Counter(),
        "mh": Counter(), "exact": defaultdict(Counter), "status": Counter(),
    } for s in samples}
    types = set()
    max_length, max_mh, rows_seen = len(model["ref"]), 0, 0
    required = {"Alias", "Type", "countEvents", "delStart", "delEnd", "delSize",
                "MH_evaluable", "junction_MH_max", "analysis_status"}
    with events_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"Missing event columns: {sorted(required - set(reader.fieldnames or []))}")
        for line, row in enumerate(reader, 2):
            rows_seen += 1
            alias = row["Alias"].strip()
            if alias not in state:
                raise ValueError(f"Events row {line}: unknown Alias {alias!r}")
            count = integer(row["countEvents"], "countEvents")
            if count < 0:
                raise ValueError(f"Events row {line}: negative countEvents")
            kind = row["Type"].strip()
            if not kind:
                raise ValueError(f"Events row {line}: missing Type")
            types.add(kind)
            if not count:
                continue
            a = state[alias]
            a["total"] += count
            a["types"][kind] += count
            flag = row["MH_evaluable"].strip().lower()
            if flag not in {"true", "false", "1", "0", ""}:
                raise ValueError(f"Events row {line}: invalid MH_evaluable")
            strict = flag in {"true", "1"}
            if kind != "DELETION":
                if strict:
                    raise ValueError(f"Events row {line}: non-DELETION marked MH-evaluable")
                continue
            a["deletion"] += count
            a["status"][row["analysis_status"]] += count
            try:
                length = integer(row["delSize"], "delSize")
                if not 0 < length <= max_length:
                    raise ValueError("delSize outside reference")
            except ValueError:
                length = None
            if length is not None:
                a["length_valid"] += count
                a["length_sum"] += count * length
                a["lengths"][length] += count
            try:
                start, end = integer(row["delStart"], "delStart"), integer(row["delEnd"], "delEnd")
                if length is None or not 0 <= start < end <= max_length or end - start != length:
                    raise ValueError("Invalid deletion coordinates")
            except ValueError:
                a["coord_invalid"] += count
                if strict:
                    raise ValueError(f"Events row {line}: invalid coordinates marked MH-evaluable")
            else:
                for position in range(start, end):
                    a["positions"][position] += count
            if strict:
                h = integer(row["junction_MH_max"], "junction_MH_max")
                if length is None or not 0 <= h <= length:
                    raise ValueError(f"Events row {line}: invalid MH length")
                a["strict"] += count
                a["strict_length_sum"] += count * length
                a["mh_sum"] += count * h
                a["mh"][h] += count
                a["exact"][length][h] += count
                max_mh = max(max_mh, h)
    if not rows_seen:
        raise ValueError("No event rows")
    tables = {k: [] for k in ["sample_metrics", "event_types_by_sample", "deletion_lengths_by_sample",
                              "deletion_positions_by_sample", "mh_distribution_by_sample", "exact_length_mh_by_sample",
                              "mh_status_by_sample"]}
    statuses = sorted({status for a in state.values() for status in a["status"]})
    for sample in samples:
        a = state[sample["Alias"]]
        total, deletion, strict = a["total"], a["deletion"], a["strict"]
        if sum(a["types"].values()) != total or sum(a["mh"].values()) != strict:
            raise AssertionError("Count conservation failed")
        metrics = dict(sample, total_reported_reads=total, pure_deletion_reads=deletion,
                       valid_length_deletion_reads=a["length_valid"],
                       invalid_length_deletion_reads=deletion-a["length_valid"],
                       invalid_coordinate_deletion_reads=a["coord_invalid"],
                       evaluable_deletion_reads=strict, unevaluable_deletion_reads=deletion-strict,
                       deletion_fraction_of_reported=fraction(deletion, total),
                       evaluable_fraction_of_deletions=fraction(strict, deletion),
                       mean_deletion_bp_all_deletions=(fraction(a["length_sum"], deletion)
                                                       if a["length_valid"] == deletion else None),
                       mean_deletion_bp_valid_lengths=fraction(a["length_sum"], a["length_valid"]),
                       mean_deletion_bp_strict_mh=fraction(a["strict_length_sum"], strict),
                       mean_mh_bp=fraction(a["mh_sum"], strict))
        for h in THRESHOLDS:
            n = sum(count for mh, count in a["mh"].items() if mh >= h)
            metrics[f"mh_ge{h}_reads"] = n
            metrics[f"mh_ge{h}_fraction_of_evaluable"] = fraction(n, strict)
        tables["sample_metrics"].append(metrics)
        for kind in sorted(types):
            tables["event_types_by_sample"].append(dict(sample, Type=kind, read_count=a["types"][kind],
                                                         total_reported_reads=total,
                                                         event_fraction=fraction(a["types"][kind], total)))
        for length in range(1, max_length + 1):
            count = a["lengths"][length]
            tables["deletion_lengths_by_sample"].append(dict(sample, deletion_length_bp=length,
                    read_count=count, pure_deletion_reads=deletion, total_reported_reads=total,
                    fraction_of_deletions=fraction(count, deletion), fraction_of_reported=fraction(count, total)))
            exact = a["exact"][length]
            n = sum(exact.values())
            record = dict(sample, deletion_length_bp=length, evaluable_reads_at_length=n,
                          evaluable_deletion_reads=strict, length_fraction_of_evaluable=fraction(n, strict))
            for h in THRESHOLDS:
                num = sum(count for mh, count in exact.items() if mh >= h)
                record[f"mh_ge{h}_reads"] = num
                record[f"mh_ge{h}_fraction"] = fraction(num, n)
            tables["exact_length_mh_by_sample"].append(record)
        for position in range(max_length):
            count = a["positions"][position]
            tables["deletion_positions_by_sample"].append(dict(sample, position_0based=position,
                    position_relative_to_junction=position-model["junction"], read_count=count,
                    pure_deletion_reads=deletion, total_reported_reads=total,
                    fraction_of_deletions=fraction(count, deletion), fraction_of_reported=fraction(count, total)))
        for h in range(max_mh + 1):
            tables["mh_distribution_by_sample"].append(dict(sample, junction_MH_bp=h,
                    read_count=a["mh"][h], evaluable_deletion_reads=strict,
                    mh_fraction=fraction(a["mh"][h], strict)))
        for status in statuses:
            tables["mh_status_by_sample"].append(dict(sample, analysis_status=status,
                    read_count=a["status"][status], pure_deletion_reads=deletion,
                    fraction_of_deletions=fraction(a["status"][status], deletion)))
    metrics = [c for c in tables["sample_metrics"][0] if c not in SAMPLE_KEYS]
    tables["group_metrics"] = group_rows(tables["sample_metrics"], [], metrics)
    groups = {
        "event_types": (["Type"], ["read_count", "event_fraction"]),
        "deletion_lengths": (["deletion_length_bp"], ["read_count", "fraction_of_deletions", "fraction_of_reported"]),
        "deletion_positions": (["position_0based", "position_relative_to_junction"], ["read_count", "fraction_of_deletions", "fraction_of_reported"]),
        "mh_distribution": (["junction_MH_bp"], ["read_count", "mh_fraction"]),
        "exact_length_mh": (["deletion_length_bp"], ["evaluable_reads_at_length", "length_fraction_of_evaluable"]
                            + [f"mh_ge{h}_fraction" for h in THRESHOLDS]),
        "mh_status": (["analysis_status"], ["read_count", "fraction_of_deletions"]),
    }
    for name, (extra, metrics) in groups.items():
        tables[name + "_by_group"] = group_rows(tables[name + "_by_sample"], extra, metrics)
    return tables, state, rows_seen


def standardize(samples, state, minimums=(20,), focus="Cy5"):
    """Descriptive direct standardization on exact-length common support.

    At each time, include lengths with >= minimum reads in *every* experiment
    in both conditions. Normalize each experiment's length counts over this
    common support, then average these distributions with equal sample weight.
    No p-values are computed; this does not adjust left/right allocation.
    """
    results, sample_results, weight_rows = [], [], []
    times = list(dict.fromkeys(s["time"] for s in samples))
    conditions = list(dict.fromkeys(s["condition"] for s in samples))
    if focus not in conditions:
        return results, sample_results, weight_rows
    for time in times:
        focus_samples = [s for s in samples if s["time"] == time and s["condition"] == focus]
        if not focus_samples:
            continue
        for other in conditions:
            if other == focus:
                continue
            comparison_samples = [s for s in samples if s["time"] == time and s["condition"] == other]
            selected = focus_samples + comparison_samples
            if not comparison_samples:
                continue
            counts = {s["Alias"]: {d: sum(c.values()) for d, c in state[s["Alias"]]["exact"].items()}
                      for s in selected}
            lengths = sorted({d for c in counts.values() for d in c})
            for minimum in minimums:
                common = [d for d in lengths if all(counts[s["Alias"]].get(d, 0) >= minimum for s in selected)]
                common_totals = {s["Alias"]: sum(counts[s["Alias"]].get(d, 0) for d in common) for s in selected}
                weights = {d: statistics.mean([counts[s["Alias"]][d] / common_totals[s["Alias"]]
                                               for s in selected]) for d in common}
                if common and not math.isclose(sum(weights.values()), 1, abs_tol=1e-12):
                    raise AssertionError("Common weights do not sum to one")
                base = dict(dataset=selected[0]["dataset"], time=time, focus=focus, comparison=other,
                            min_reads_per_length=minimum, n_common_lengths=len(common))
                for d, weight in weights.items():
                    weight_rows.append(dict(base, deletion_length_bp=d, common_weight_fraction=weight))
                for threshold in THRESHOLDS:
                    local = []
                    for sample in selected:
                        alias = sample["Alias"]
                        adjusted = sum(weights[d] * sum(n for h, n in state[alias]["exact"][d].items() if h >= threshold)
                                       / counts[alias][d] for d in common) if common else None
                        coverage = fraction(common_totals[alias], state[alias]["strict"])
                        local.append(dict(base, mh_threshold=threshold, Alias=alias, condition=sample["condition"],
                                          replicate=sample["replicate"], adjusted_fraction=adjusted,
                                          coverage_fraction=coverage, common_support_reads=common_totals[alias],
                                          evaluable_deletion_reads=state[alias]["strict"]))
                    sample_results.extend(local)
                    fg = [r for r in local if r["condition"] == focus]
                    cg = [r for r in local if r["condition"] == other]
                    fm, fs, fn = summary([r["adjusted_fraction"] for r in fg])
                    cm, cs, cn = summary([r["adjusted_fraction"] for r in cg])
                    results.append(dict(base, mh_threshold=threshold, common_lengths=";".join(map(str, common)),
                                        focus_adjusted_fraction_mean=fm, focus_adjusted_fraction_sd=fs, focus_n_valid=fn,
                                        comparison_adjusted_fraction_mean=cm, comparison_adjusted_fraction_sd=cs,
                                        comparison_n_valid=cn, difference_fraction=fm-cm if fm is not None and cm is not None else None,
                                        focus_coverage_fraction_mean=summary([r["coverage_fraction"] for r in fg])[0],
                                        comparison_coverage_fraction_mean=summary([r["coverage_fraction"] for r in cg])[0],
                                        status="ok" if common else "no_common_lengths"))
    return results, sample_results, weight_rows


def run(events_path, samples_path, dataset, config_path, outdir, minimums=(20,)):
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise ValueError("Output folder must be new or empty")
    if any(integer(m, "min_reads") <= 0 for m in minimums):
        raise ValueError("min_reads must be positive integers")
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    model = prepare_model(config)
    samples = read_samples(samples_path, dataset)
    tables, state, rows_seen = summarize(events_path, samples, model)
    if dataset == "modify":
        result, sample_result, weights = standardize(samples, state, sorted(set(minimums)))
        tables.update(length_standardized=result, length_standardized_samples=sample_result,
                      length_standardization_weights=weights)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, records in tables.items():
        if records:
            write_csv(outdir / (name + ".csv"), records)
    metadata = {
        "script_version": VERSION, "dataset": dataset, "n_samples": len(samples), "event_rows": rows_seen,
        "input_sha256": hashlib.sha256(events_path.read_bytes()).hexdigest(),
        "samples_sha256": hashlib.sha256(samples_path.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "total_reported_reads": sum(a["total"] for a in state.values()),
        "pure_deletion_reads": sum(a["deletion"] for a in state.values()),
        "evaluable_deletion_reads": sum(a["strict"] for a in state.values()),
        "unevaluable_deletion_reads": sum(a["deletion"]-a["strict"] for a in state.values()),
        "rate_unit": "fraction, not percent", "group_unit": "independent experiment; equal-weight mean and sample SD (ddof=1)",
        "missing": "Blank numeric cells mean unavailable; zero denominators yield NA, not zero.",
        "all_deletion_population": "Positive countEvents with Type exactly DELETION, including internal-only and Raw-mismatch events.",
        "strict_mh_population": "Positive-count DELETION events with MH_evaluable=True from analyze_junction_mh.py.",
        "positions": "0-based reference base positions; membership delStart <= position < delEnd, using original SIQ coordinates. Equivalent alignments are not redistributed.",
        "all_deletion_length_mean": "Count-weighted per sample. Unavailable if any DELETION length is invalid; a separate valid-length mean is also reported.",
        "standardization": "Within each time and Cy5 comparison, common lengths meet minimum in every sample. Common weights are equally averaged normalized sample length frequencies. Descriptive only; left/right allocation is not matched.",
        "minimum_reads": list(minimums),
        "tables": {name: len(records) for name, records in tables.items()},
    }
    (outdir / "summary_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Summarized {rows_seen:,} events for {len(samples)} {dataset} samples: {outdir}")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--dataset", required=True, choices=["modify", "inhibitor"])
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-reads", nargs="+", type=int, default=[20], help="Per-sample common-length minimum(s); default 20")
    args = parser.parse_args()
    try:
        run(args.events, args.samples, args.dataset, args.config, args.outdir, args.min_reads)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
