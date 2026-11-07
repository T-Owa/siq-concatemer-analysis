#!/usr/bin/env python3
"""Extract minimal, portable SIQ input tables without changing the workbook.

Only the fields required by the downstream analysis are exported. Original
file paths and other unrelated workbook columns are not copied. Requires
openpyxl; subsequent analysis of the CSV/GZ files uses the standard library.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

from analyze_junction_mh import (
    digest, dna, integer, prepare_model, validate_runinfo, worksheet_records,
)

RAW_FIELDS = [
    "source_excel_row", "Alias", "Type", "countEvents", "Raw",
    "delStart", "delEnd", "delSize", "homology", "homologyLength",
]
SAMPLE_FIELDS = ["dataset", "Alias", "condition", "condition_label", "time", "replicate"]
LABELS = {
    "ds": "Unmodified", "phos": "Phos", "Bio": "Biotin", "dsp": "dSp",
    "sp18": "Sp18", "Bractppin": "Bractoppin", "DI03": "D-I03",
}


def parse_alias(alias, dataset):
    try:
        if dataset == "modify":
            condition, time, replicate = alias.rsplit("-", 2)
        elif dataset == "inhibitor":
            time, remainder = alias.split("-", 1)
            condition, replicate = remainder.rsplit("-", 1)
        else:
            raise ValueError("dataset must be modify or inhibitor")
        replicate = integer(replicate, "replicate")
        if not condition or time not in {"12h", "24h"} or replicate < 1:
            raise ValueError("invalid condition, time or replicate")
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Cannot parse {dataset} Alias {alias!r}: {exc}") from None
    return dict(dataset=dataset, Alias=alias, condition=condition,
                condition_label=LABELS.get(condition, condition),
                time=time, replicate=replicate)


def ensure_empty(path):
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"Output folder is not empty: {path}")


def write_table(path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def extract_workbook(input_path, dataset, config_path, destination):
    """Write one dataset into an existing empty staging directory."""
    import openpyxl

    input_path, config_path, destination = map(Path, (input_path, config_path, destination))
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    model = prepare_model(config)
    warnings.filterwarnings("ignore", message="Workbook contains no default style")
    workbook = openpyxl.load_workbook(input_path, read_only=True, data_only=True)
    try:
        checked = validate_runinfo(workbook, config, model)
        samples = [parse_alias(r["Alias"], dataset) for r in checked]
        aliases = {r["Alias"] for r in samples}
        if "rawData" not in workbook.sheetnames:
            raise ValueError("No rawData sheet in input workbook")
        output_path = destination / f"{dataset}_events.csv.gz"
        rows_seen = 0
        total_reads = 0
        source_required = set(RAW_FIELDS) - {"source_excel_row", "homology", "homologyLength"}
        # Explicit filename='' and mtime=0 make gzip bytes deterministic.
        with output_path.open("wb") as raw_handle:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8-sig", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=RAW_FIELDS)
                    writer.writeheader()
                    for row_number, row in worksheet_records(workbook["rawData"]):
                        if not source_required <= row.keys():
                            raise ValueError(f"Missing required columns: {sorted(source_required - row.keys())}")
                        alias = str(row["Alias"] or "").strip()
                        if alias not in aliases:
                            raise ValueError(f"rawData row {row_number}: unknown Alias {alias!r}")
                        count = integer(row["countEvents"], "countEvents")
                        if count < 0:
                            raise ValueError(f"rawData row {row_number}: negative countEvents")
                        kind = str(row["Type"] or "").strip()
                        raw = dna(row["Raw"])
                        if kind == "WT" and raw != model["ref"]:
                            raise ValueError(f"rawData row {row_number}: WT Raw differs from reference")
                        record = {k: row.get(k) for k in RAW_FIELDS}
                        record.update(source_excel_row=row_number, Alias=alias,
                                      Type=kind, countEvents=count, Raw=raw)
                        writer.writerow(record)
                        rows_seen += 1
                        total_reads += count
        if not rows_seen:
            raise ValueError("No rawData rows")
        settings = []
        if "SIQRun" in workbook.sheetnames:
            # Explicit allowlist avoids exporting unrelated columns or paths.
            settings = [{"Info": r.get("Info"), "Value": r.get("Value")}
                        for _, r in worksheet_records(workbook["SIQRun"])]
        metadata = {
            "dataset": dataset,
            "source_filename": input_path.name,
            "source_sha256": digest(input_path),
            "rawData_rows": rows_seen,
            "n_samples": len(samples),
            "total_reported_reads": total_reads,
            "export_filename": output_path.name,
            "export_sha256": digest(output_path),
            "export_fields": RAW_FIELDS,
            "SIQRun": settings,
            "runinfo_checked": checked,
            "config": config,
            "config_sha256": digest(config_path),
            "export_policy": "All rawData rows, including zero-count rows. Only the listed required analysis columns are retained. Original file paths are omitted.",
        }
        return samples, metadata
    finally:
        workbook.close()


def run(input_path, dataset, config_path, outdir):
    ensure_empty(outdir)
    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="siq_import_", dir=outdir.parent) as temp:
        staging = Path(temp)
        samples, metadata = extract_workbook(input_path, dataset, config_path, staging)
        write_table(staging / "samples.csv", SAMPLE_FIELDS, samples)
        (staging / "source_metadata.json").write_text(
            json.dumps({"datasets": [metadata]}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        ensure_empty(outdir)
        outdir.mkdir(exist_ok=True)
        for path in staging.iterdir():
            shutil.move(str(path), str(outdir / path.name))
    print(f"Imported {metadata['rawData_rows']:,} rows for {len(samples)} samples into {outdir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--dataset", required=True, choices=["modify", "inhibitor"])
    parser.add_argument("--config", type=Path, default=Path(__file__).resolve().parents[1] / "config/junction.json")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        run(args.input.resolve(), args.dataset, args.config.resolve(), args.outdir.resolve())
    except (ValueError, OSError, ImportError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
