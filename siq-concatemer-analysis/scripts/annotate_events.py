#!/usr/bin/env python3
"""Apply the original strict junction-MH calculation to portable SIQ CSV input.

Scientific annotation is delegated unchanged to analyze_junction_mh.py.
The original observed Raw sequence is required to agree with the SIQ-coordinate
single-deletion reconstruction. Unevaluable records retain blank MH values.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

from analyze_junction_mh import (
    EVENT_FIELDS, VERSION, annotate_deletion, digest, dna, integer, prepare_model,
)


def open_text(path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8-sig", newline="")


def ensure_empty(path):
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"Output folder is not empty: {path}")


def read_samples(path, dataset):
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        required = {"dataset", "Alias", "condition", "condition_label", "time", "replicate"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"Missing sample columns: {sorted(required - set(reader.fieldnames or []))}")
        samples = [r for r in reader if r["dataset"] == dataset]
    aliases = [r["Alias"].strip() for r in samples]
    if not aliases or any(not a for a in aliases) or len(aliases) != len(set(aliases)):
        raise ValueError("Selected sample aliases must be nonempty and unique")
    for row in samples:
        row["Alias"] = row["Alias"].strip()
        rep = integer(row["replicate"], "replicate")
        if rep < 1 or not row["condition"] or not row["time"]:
            raise ValueError(f"Invalid sample metadata for {row['Alias']!r}")
    return samples


def run(input_path, samples_path, dataset, config_path, outdir):
    ensure_empty(outdir)
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    model = prepare_model(config)
    samples = read_samples(samples_path, dataset)
    aliases = {r["Alias"] for r in samples}
    outdir.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    status_reads = Counter()
    status_rows = Counter()
    observed_aliases = set()
    source_rows = set()
    with tempfile.TemporaryDirectory(prefix="siq_annotate_", dir=outdir.parent) as temp:
        staging = Path(temp)
        with open_text(input_path) as source, (staging / "events.csv").open("w", encoding="utf-8-sig", newline="") as target:
            reader = csv.DictReader(source)
            required = {"source_excel_row", "Alias", "Type", "countEvents", "Raw", "delStart", "delEnd", "delSize"}
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f"Missing input columns: {sorted(required - set(reader.fieldnames or []))}")
            writer = csv.DictWriter(target, fieldnames=EVENT_FIELDS, extrasaction="ignore")
            writer.writeheader()
            for csv_row, row in enumerate(reader, start=2):
                rownum = integer(row["source_excel_row"], "source_excel_row")
                if rownum < 2 or rownum in source_rows:
                    raise ValueError(f"CSV row {csv_row}: invalid or duplicate source_excel_row {rownum}")
                source_rows.add(rownum)
                alias = str(row["Alias"] or "").strip()
                if alias not in aliases:
                    raise ValueError(f"CSV row {csv_row}: Alias {alias!r} missing from selected manifest")
                observed_aliases.add(alias)
                count = integer(row["countEvents"], "countEvents")
                if count < 0:
                    raise ValueError(f"CSV row {csv_row}: negative countEvents")
                kind, raw = str(row["Type"] or "").strip(), dna(row["Raw"])
                if kind == "WT" and raw != model["ref"]:
                    raise ValueError(f"CSV row {csv_row}: WT Raw differs from configured reference")
                event = {k: row.get(k) for k in ("delStart", "delEnd", "delSize")}
                event.update(source_excel_row=rownum, Alias=alias, Type=kind,
                             countEvents=count, Raw=raw, siq_homology=row.get("homology"),
                             siq_homologyLength=row.get("homologyLength"), MH_evaluable=False)
                counts["total_reported_reads"] += count
                if count == 0:
                    event["analysis_status"] = "zero_count_reference"
                elif kind != "DELETION":
                    event["analysis_status"] = "not_target_type"
                else:
                    counts["pure_deletion_reads"] += count
                    event.update(annotate_deletion(raw, row["delStart"], row["delEnd"], row["delSize"], model))
                    if event["MH_evaluable"]:
                        counts["evaluable_deletion_reads"] += count
                    else:
                        counts["unevaluable_deletion_reads"] += count
                status_reads[event["analysis_status"]] += count
                status_rows[event["analysis_status"]] += 1
                counts["rawData_rows"] += 1
                writer.writerow(event)
        if not counts["rawData_rows"]:
            raise ValueError("No input events")
        if observed_aliases != aliases:
            raise ValueError(f"Manifest samples without rows: {sorted(aliases - observed_aliases)}")
        metadata = {
            "dataset": dataset, "core_version": VERSION,
            "python_version": platform.python_version(),
            "input_filename": input_path.name, "input_sha256": digest(input_path),
            "samples_filename": samples_path.name, "samples_sha256": digest(samples_path),
            "config_sha256": digest(config_path), "config": config,
            "core_sha256": digest(Path(__file__).with_name("analyze_junction_mh.py")),
            "n_samples": len(samples), **dict(counts),
            "status_reads": dict(sorted(status_reads.items())),
            "status_rows": dict(sorted(status_rows.items())),
            "definition": "Maximum junction-compatible homology from the complete observed Raw under the exact single-deletion model; not proof of biochemical MH use.",
            "denominator": "Strict MH summaries use positive-count pure DELETION records with MH_evaluable=True; blank MH is not MH=0.",
        }
        (staging / "annotation_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        ensure_empty(outdir)
        outdir.mkdir(exist_ok=True)
        for path in staging.iterdir():
            shutil.move(str(path), str(outdir / path.name))
    print(f"Annotated {counts['rawData_rows']:,} rows for {len(samples)} samples; {counts['evaluable_deletion_reads']:,} evaluable deletion reads")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--dataset", required=True, choices=["modify", "inhibitor"])
    parser.add_argument("--config", type=Path, default=Path(__file__).resolve().parents[1] / "config/junction.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        run(args.input.resolve(), args.samples.resolve(), args.dataset, args.config.resolve(), args.outdir.resolve())
    except (ValueError, OSError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
