#!/usr/bin/env python3
"""Aggregate coordinate-reconstructed SIQ deletion alleles and per-group Top5.

All-event denominators are the sum of countEvents in the supplied SIQ event
export. They are not FASTQ, merged-read, or independently measured molecule
counts. Positive pure DELETION reads remain in denominators even if Raw does
not reconstruct exactly or no junction-compatible placement exists.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from analyze_junction_mh import annotate_deletion, prepare_model, integer, dna

VERSION = "1.0.0"


def load_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def hash_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mean_sd(values):
    valid = [x for x in values if x is not None and math.isfinite(x)]
    return (statistics.mean(valid) if valid else None,
            statistics.stdev(valid) if len(valid) > 1 else None, len(valid))


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def write_csv(path, records, empty_fields):
    fields = list(records[0]) if records else empty_fields
    with Path(path).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in records:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict))
                             else v for k, v in row.items()})


def sample_metadata(rows, dataset):
    samples = {}
    for row in rows:
        if row.get("dataset") != dataset:
            continue
        for field in ("Alias", "condition", "condition_label", "time", "replicate"):
            if not row.get(field):
                raise ValueError(f"Missing sample metadata: {field}")
        alias = row["Alias"]
        if alias in samples:
            raise ValueError(f"Duplicate sample Alias: {alias}")
        samples[alias] = {k: row[k] for k in
                          ("dataset", "Alias", "condition", "condition_label", "time")}
        samples[alias]["replicate"] = integer(row["replicate"], "replicate")
    if not samples:
        raise ValueError(f"No sample metadata for dataset {dataset!r}")
    group_replicates = set()
    for s in samples.values():
        key = (s["condition"], s["time"], s["replicate"])
        if key in group_replicates:
            raise ValueError(f"Duplicate condition/time/replicate: {key}")
        group_replicates.add(key)
    return samples


def catalog_by_allele(rows, ref):
    catalog, identifiers = {}, set()
    for row in rows:
        allele = dna(row["reconstructed_allele"])
        pid = row["pattern_id"]
        if allele in catalog or pid in identifiers:
            raise ValueError("Catalog contains duplicate allele or pattern ID")
        s = integer(row["canonical_start"], "catalog canonical_start")
        e = integer(row["canonical_end"], "catalog canonical_end")
        d = integer(row["length"], "catalog length")
        if not (0 <= s < e <= len(ref) and e - s == d and ref[:s] + ref[e:] == allele):
            raise ValueError(f"Catalog {pid}: allele and coordinates disagree")
        catalog[allele] = row
        identifiers.add(pid)
    return catalog


def interval_distance(s, d, junction):
    return max(s - junction, junction - (s + d), 0)


def signed(value):
    return f"+{value}" if value > 0 else str(value)


def make_definition(allele, coordinates, catalog, model, total, mismatch):
    ref, j = model["ref"], model["junction"]
    d = len(ref) - len(allele)
    canonical_start = max(s for s, _ in coordinates)
    canonical_end = canonical_start + d
    known = catalog.get(allele, {})
    pid = known.get("pattern_id") or "ALT-" + hashlib.sha256(allele.encode("ascii")).hexdigest()[:16]
    starts = [s for s in range(len(ref) - d + 1) if ref[:s] + ref[s + d:] == allele]
    display = min(starts, key=lambda s: (interval_distance(s, d, j),
                                       abs(s - canonical_start), s))
    annotated = annotate_deletion(allele, canonical_start, canonical_end, d, model)
    valid = [s for s in starts if model["left_primer_end"] <= s <= j <= s + d <= model["right_primer_start"]]
    return {
        "pattern_id": pid, "legacy_id": known.get("legacy_id", ""), "length": d,
        "canonical_start": canonical_start, "canonical_end": canonical_end,
        "display_start": display, "display_end": display + d,
        "display_start_relative": display - j, "display_end_relative": display + d - j,
        "display_interval": f"[{signed(display-j)},{signed(display+d-j)})",
        "source_coordinate_pairs": [list(x) for x in sorted(coordinates)],
        "equivalent_starts": starts, "junction_starts": valid,
        "junction_compatible": bool(valid),
        "internal_alternative": any(not s <= j <= s + d for s in starts),
        "classification": annotated["analysis_status"],
        "mh_length": annotated.get("junction_MH_max"),
        "mh_motif": annotated.get("junction_MH_sequence", ""),
        "mh_annotation_method": "SIQ reconstructed reference allele; maximum junction-compatible MH",
        "tail_loss_min": annotated.get("tail_loss_min"),
        "tail_loss_max": annotated.get("tail_loss_max"),
        "head_loss_min": annotated.get("head_loss_min"),
        "head_loss_max": annotated.get("head_loss_max"),
        "tail_head_pairs": [[j-s, s+d-j] for s in valid],
        "deleted_sequence": ref[display:display+d],
        "reconstructed_allele": allele, "total_count": total,
        "raw_mismatch_count": mismatch,
    }


def aggregate(events, sample_rows, catalog_rows, config, dataset):
    """Return all output tables without writing files; all fractions are 0..1."""
    model = prepare_model(config)
    ref = model["ref"]
    samples = sample_metadata(sample_rows, dataset)
    catalog = catalog_by_allele(catalog_rows, ref)
    counts, mismatches, strict_counts = Counter(), Counter(), Counter()
    totals = {a: Counter() for a in samples}
    length_totals = {(a, d): Counter() for a in samples for d in (6, 12)}
    coords = defaultdict(set)
    allele_totals, allele_mismatches = Counter(), Counter()
    annotations = {}
    required = {"Alias", "Type", "countEvents", "Raw", "delStart", "delEnd", "delSize"}
    n_rows = 0
    for row in events:
        n_rows += 1
        if not required <= row.keys():
            raise ValueError(f"Missing events columns: {sorted(required - row.keys())}")
        alias = row["Alias"].strip()
        if alias not in samples:
            raise ValueError(f"Event Alias not in selected sample metadata: {alias}")
        n = integer(row["countEvents"], "countEvents")
        if n < 0:
            raise ValueError("Negative countEvents")
        totals[alias]["all_event_reads"] += n
        if n == 0 or row["Type"].strip() != "DELETION":
            continue
        totals[alias]["deletion_reads"] += n
        try:
            d = integer(row["delSize"], "delSize")
        except ValueError:
            d = None
        buckets = [totals[alias]]
        if d in (6, 12):
            buckets.append(length_totals[alias, d])
            length_totals[alias, d]["deletion_reads"] += n
        try:
            s, e = integer(row["delStart"], "delStart"), integer(row["delEnd"], "delEnd")
            good = 0 <= s < e <= len(ref) and e - s == d
        except ValueError:
            good = False
        if not good:
            for bucket in buckets:
                bucket["unclassified_count"] += n
            continue
        allele = ref[:s] + ref[e:]
        coords[allele].add((s, e))
        counts[alias, allele] += n
        allele_totals[allele] += n
        raw = dna(row["Raw"])
        mismatch = raw != allele
        key = (raw, s, e, d)
        if key not in annotations:
            annotations[key] = annotate_deletion(raw, s, e, d, model)
        strict = annotations[key]["MH_evaluable"]
        reconstructed_key = (allele, s, e, d)
        if reconstructed_key not in annotations:
            annotations[reconstructed_key] = annotate_deletion(allele, s, e, d, model)
        reconstructed = annotations[reconstructed_key]
        internal_only = reconstructed["junction_candidate_count"] == 0
        if mismatch:
            mismatches[alias, allele] += n
            allele_mismatches[allele] += n
        if strict:
            strict_counts[alias, allele] += n
        for bucket in buckets:
            bucket["classified_count"] += n
            bucket["raw_mismatch_count"] += n * mismatch
            bucket["strict_evaluable_denominator"] += n * strict
            bucket["internal_only_count"] += n * internal_only
    if not n_rows:
        raise ValueError("Event CSV contains no records")
    definitions = [make_definition(a, coords[a], catalog, model, allele_totals[a], allele_mismatches[a])
                   for a in coords]
    definitions.sort(key=lambda x: (x["length"], x["pattern_id"]))
    if len({x["pattern_id"] for x in definitions}) != len(definitions):
        raise ValueError("Pattern ID collision")
    by_allele = {r["reconstructed_allele"]: r for r in definitions}
    groups = defaultdict(list)
    for alias, sample in samples.items():
        groups[sample["condition"], sample["time"]].append(alias)
    for aliases in groups.values():
        aliases.sort(key=lambda a: (samples[a]["replicate"], a))
    tables = {"pattern_definitions": definitions}
    scopes = [("all", None)] + ([("length", 6), ("length", 12)] if dataset == "modify" else [])
    for prefix, length in scopes:
        relevant = {a for a in coords if length is None or len(ref)-len(a) == length}
        sample_totals_rows, sample_rows_out, group_rows, top_rows, coverage_rows = [], [], [], [], []
        for alias, meta in samples.items():
            total = totals[alias] if length is None else length_totals[alias, length]
            assert total["deletion_reads"] == total["classified_count"] + total["unclassified_count"]
            sample_totals_rows.append({
                "dataset": dataset, "length": length, **{k: meta[k] for k in
                ("condition", "condition_label", "time", "replicate")}, "alias": alias,
                "denominator": total["deletion_reads"],
                "all_deletion_denominator": totals[alias]["deletion_reads"],
                "all_event_denominator": totals[alias]["all_event_reads"],
                "strict_evaluable_denominator": total["strict_evaluable_denominator"],
                "raw_mismatch_count": total["raw_mismatch_count"],
                "internal_only_count": total["internal_only_count"],
                "unclassified_count": total["unclassified_count"],
            })
        total_lookup = {r["alias"]: r for r in sample_totals_rows}
        for (condition, time), aliases in sorted(groups.items()):
            # Exact-length output includes every observed pattern in every group.
            # All-length output omits patterns absent from the entire group.
            alleles = relevant if length is not None else {
                a for a in relevant if any(counts[alias, a] > 0 for alias in aliases)}
            local_rows = []
            for allele in sorted(alleles, key=lambda a: by_allele[a]["pattern_id"]):
                definition = by_allele[allele]
                reps = []
                for alias in aliases:
                    t = total_lookup[alias]
                    n = counts[alias, allele]
                    record = {**t, "pattern_id": definition["pattern_id"],
                              "legacy_id": definition["legacy_id"],
                              "pattern_length": definition["length"], "count": n,
                              "fraction": ratio(n, t["denominator"]),
                              "fraction_all_events": ratio(n, t["all_event_denominator"]),
                              "fraction_all_deletions": ratio(n, t["all_deletion_denominator"]),
                              "pattern_raw_mismatch_count": mismatches[alias, allele],
                              "strict_evaluable_count": strict_counts[alias, allele],
                              "strict_evaluable_fraction": ratio(strict_counts[alias, allele], t["strict_evaluable_denominator"])}
                    sample_rows_out.append(record)
                    reps.append(record)
                row = {"dataset": dataset, "length": length, "condition": condition,
                       "condition_label": samples[aliases[0]]["condition_label"], "time": time,
                       "pattern_id": definition["pattern_id"], "legacy_id": definition["legacy_id"],
                       "pattern_length": definition["length"], "mh_length": definition["mh_length"],
                       "junction_compatible": definition["junction_compatible"], "n_samples": len(aliases)}
                for field in ("fraction", "fraction_all_events", "fraction_all_deletions", "strict_evaluable_fraction"):
                    m, sd, n_valid = mean_sd([r[field] for r in reps])
                    row["mean_" + field], row["sd_" + field], row["n_valid_" + field] = m, sd, n_valid
                row["n_valid"] = row["n_valid_fraction"]
                row["n"] = row["n_valid"]
                row.update(replicate_aliases=aliases, replicate_counts=[r["count"] for r in reps],
                           denominators=[r["denominator"] for r in reps],
                           all_event_denominators=[r["all_event_denominator"] for r in reps],
                           all_deletion_denominators=[r["all_deletion_denominator"] for r in reps],
                           replicate_fractions=[r["fraction"] for r in reps],
                           raw_mismatch_counts=[r["pattern_raw_mismatch_count"] for r in reps],
                           strict_evaluable_counts=[r["strict_evaluable_count"] for r in reps],
                           strict_evaluable_denominators=[r["strict_evaluable_denominator"] for r in reps],
                           strict_evaluable_fractions=[r["strict_evaluable_fraction"] for r in reps])
                local_rows.append(row)
            local_rows.sort(key=lambda r: (-(r["mean_fraction"] if r["mean_fraction"] is not None else -1), r["pattern_id"]))
            positive_rank = 0
            for row in local_rows:
                positive = row["mean_fraction"] is not None and row["mean_fraction"] > 0
                if positive:
                    positive_rank += 1
                row["rank"] = positive_rank if positive else None
                row["top5"] = positive and positive_rank <= 5
                if row["top5"]:
                    top_rows.append(row)
            group_rows.extend(local_rows)
            selected = [r for r in local_rows if r["top5"]]
            top_counts = [sum(r["replicate_counts"][i] for r in selected) for i in range(len(aliases))]
            denoms = [total_lookup[a]["denominator"] for a in aliases]
            fractions = [ratio(n, d) for n, d in zip(top_counts, denoms)]
            other_counts = [d-n for n, d in zip(top_counts, denoms)]
            other_fracs = [ratio(n, d) for n, d in zip(other_counts, denoms)]
            m, sd, n_valid = mean_sd(fractions)
            om, osd, _ = mean_sd(other_fracs)
            coverage_rows.append({
                "dataset": dataset, "length": length, "condition": condition,
                "condition_label": samples[aliases[0]]["condition_label"], "time": time,
                "n_top_patterns": len(selected), "top_pattern_ids": [r["pattern_id"] for r in selected],
                "replicate_aliases": aliases, "top5_counts": top_counts,
                "denominators": denoms, "top5_fractions": fractions,
                "top5_mean_fraction": m, "top5_sd_fraction": sd, "n_valid": n_valid,
                "other_mean_fraction": om, "other_sd_fraction": osd,
                "other_counts": other_counts, "other_fractions": other_fracs,
            })
        for suffix, rows in (("sample_totals", sample_totals_rows), ("sample_patterns", sample_rows_out),
                             ("group_patterns", group_rows), ("top5", top_rows), ("coverage", coverage_rows)):
            tables.setdefault(prefix + "_" + suffix, []).extend(rows)
    return tables


def run(events_path, samples_path, dataset, config_path, catalog_path, outdir):
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise ValueError("Output folder must be new or empty")
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    tables = aggregate(load_csv(events_path), load_csv(samples_path), load_csv(catalog_path), config, dataset)
    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="patterns_", dir=outdir.parent) as tmp:
        stage = Path(tmp)
        for name, rows in tables.items():
            write_csv(stage / (name + ".csv"), rows, ["dataset", "length", "pattern_id"])
        metadata = {"script_version": VERSION, "dataset": dataset,
                    "input_sha256": {"events": hash_file(events_path), "samples": hash_file(samples_path),
                                     "config": hash_file(config_path), "catalog": hash_file(catalog_path)},
                    "row_counts": {k: len(v) for k, v in tables.items()},
                    "fractions": "0..1; mean of per-sample fractions; sample SD (ddof=1)",
                    "ranking": "positive group mean fraction within selected pure-DELETION denominator, descending; pattern_id tie-break",
                    "mh_definition": "Coordinate-reconstructed allele annotation, not necessarily strict full-Raw-evaluable MH",
                    "all_event_denominator": "sum of all countEvents in SIQ output, not MergedCorrect or original sequencing reads",
                    "display_coordinates": "Equivalent placement nearest junction, then minimum shift from maximum source delStart; actual cleavage positions not inferred",
                    "missing_denominator": "blank/NA; not zero; n_valid records available samples"}
        (stage / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        if outdir.exists() and any(outdir.iterdir()):
            raise ValueError("Output folder became nonempty; not overwritten")
        outdir.mkdir(exist_ok=True)
        for path in stage.iterdir():
            shutil.move(str(path), str(outdir / path.name))
    print(f"Patterns: {dataset}; {len(tables['pattern_definitions'])} reconstructed alleles; output {outdir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--dataset", choices=["modify", "inhibitor"], required=True)
    parser.add_argument("--config", type=Path, default=Path(__file__).resolve().parents[1] / "config/junction.json")
    parser.add_argument("--catalog", type=Path, default=Path(__file__).resolve().parents[1] / "config/pattern_catalog.csv")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        run(args.events, args.samples, args.dataset, args.config, args.catalog, args.outdir)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
