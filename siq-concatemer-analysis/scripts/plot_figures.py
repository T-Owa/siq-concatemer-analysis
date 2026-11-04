#!/usr/bin/env python3
"""Create optional static publication figures from the generated CSV tables.

These figures reproduce the numerical analysis, not the historical Figure
Studio interface. MH means maximum sequence-compatible junction homology.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

MOD_ORDER = ["ds", "phos", "Acryd", "AzideN", "AmMC6", "Cy5", "Bio", "InvddT", "dsp", "sp18"]
TYPE_ORDER = ["WT", "DELETION", "INSERTION", "SNV", "DELINS", "TINS", "DELINS_SNV",
              "TANDEMDUPLICATION", "TANDEMDUPLICATION_COMPOUND"]


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    if value is None or str(value).strip() in ("", "None", "NA", "nan"):
        return float("nan")
    return float(value)


def truth(value):
    return str(value).lower() in ("true", "1")


def conditions_for(rows, requested=None):
    present = list(dict.fromkeys(r["condition"] for r in rows))
    if requested:
        missing = set(requested) - set(present)
        if missing:
            raise ValueError(f"Requested conditions are absent: {sorted(missing)}")
        return requested
    return [x for x in MOD_ORDER if x in present] + sorted(set(present) - set(MOD_ORDER))


def label_map(rows):
    return {r["condition"]: r.get("condition_label", r["condition"]) for r in rows}


def save(fig, stem, args):
    for extension in args.formats:
        buffer = io.BytesIO()
        fig.savefig(buffer, format=extension, dpi=args.dpi,
                    facecolor="white", bbox_inches="tight", pad_inches=0.08)
        payload = buffer.getvalue()
        if not payload:
            raise RuntimeError(f"Figure export was empty: {stem}.{extension}")
        (args.outdir / f"{stem}.{extension}").write_bytes(payload)
    plt.close(fig)


def plot_top5(pattern_dir, args):
    definitions = {r["pattern_id"]: r for r in read_rows(pattern_dir / "pattern_definitions.csv")}
    rows = read_rows(pattern_dir / "length_group_patterns.csv")
    coverage = read_rows(pattern_dir / "length_coverage.csv")
    names = label_map(rows)
    conditions = conditions_for(rows, args.conditions)
    created = []
    for length in args.lengths:
        for time in args.times:
            subset = [r for r in rows if int(r["length"]) == length and r["time"] == time
                      and r["condition"] in conditions]
            if not subset:
                continue
            ids = sorted({r["pattern_id"] for r in subset if truth(r["top5"])})
            if not ids:
                continue
            lookup = {(r["pattern_id"], r["condition"]): r for r in subset}
            matrix = np.array([[100 * number(lookup.get((pid, c), {}).get("mean_fraction"))
                                for c in conditions] for pid in ids])
            starts = [number(definitions[pid]["display_start_relative"]) for pid in ids]
            ends = [number(definitions[pid]["display_end_relative"]) for pid in ids]
            limit = max(10, 5 * math.ceil(max([abs(x) for x in starts + ends]) / 5) + 5)
            fig = plt.figure(figsize=(args.width_mm / 25.4, max(4.0, len(ids) * 0.38 + 1.9)))
            grid = fig.add_gridspec(1, 3, width_ratios=(2.5, 1.45, max(4.2, len(conditions) * 0.48)),
                                   wspace=0.13, left=0.11, right=0.96, bottom=0.32, top=0.83)
            ax_map, ax_meta, ax_heat = [fig.add_subplot(grid[0, i]) for i in range(3)]
            y = np.arange(len(ids))
            ax_map.axvspan(-limit, 0, color="#dce9f2", zorder=0)
            ax_map.axvspan(0, limit, color="#f8eadb", zorder=0)
            ax_map.axvline(0, color="#777777", ls="--", lw=0.8)
            row_labels = []
            for index, pid in enumerate(ids):
                d = definitions[pid]
                internal = not truth(d["junction_compatible"])
                row_labels.append(pid + ("†" if internal else ""))
                ax_map.plot([-limit, limit], [index, index], color="#c0c6cb", lw=0.8, zorder=1)
                ax_map.add_patch(Rectangle((starts[index], index - 0.20), ends[index] - starts[index],
                                           0.40, facecolor="#bf5b35" if internal else "#314e68",
                                           edgecolor="none", zorder=2))
                h = number(d.get("mh_length"))
                mh = "MH: NA" if not math.isfinite(h) else f"MH: {int(h)} bp"
                motif = d.get("mh_motif", "")
                if motif:
                    mh += " " + motif
                interval = d.get("display_interval", f"[{starts[index]:g},{ends[index]:g})")
                ax_meta.text(0, index - 0.20, interval, va="center", fontsize=7.2)
                ax_meta.text(0, index + 0.20, mh, va="center", fontsize=7.2,
                             color="#777777" if internal else "#333333")
            for ax in (ax_map, ax_meta, ax_heat):
                ax.set_ylim(len(ids) - 0.5, -0.5)
            ax_map.set_xlim(-limit, limit)
            ax_map.set_yticks(y, row_labels, fontsize=8)
            ax_map.tick_params(axis="y", length=0)
            ax_map.set_xticks([-limit, 0, limit], [str(-limit), "0", f"+{limit}"])
            ax_map.set_xlabel("Position relative to junction", fontsize=8)
            ax_map.set_title("Tail           Head", fontsize=9, pad=14)
            ax_meta.set_xlim(0, 1)
            ax_meta.axis("off")
            ax_meta.set_title("Interval / MH candidate", fontsize=8, loc="left", pad=14)
            palette = plt.get_cmap(args.cmap).copy()
            palette.set_bad("#ececec")
            im = ax_heat.imshow(np.ma.masked_invalid(matrix), vmin=0, vmax=100, cmap=palette,
                                aspect="auto", interpolation="none")
            ax_heat.set_yticks([])
            ax_heat.set_xticks(np.arange(len(conditions)), [names[x] for x in conditions],
                              rotation=55, ha="left", rotation_mode="anchor", fontsize=8)
            ax_heat.xaxis.tick_top()
            ax_heat.tick_params(axis="both", length=0)
            ax_heat.set_xticks(np.arange(len(conditions) + 1) - 0.5, minor=True)
            ax_heat.set_yticks(np.arange(len(ids) + 1) - 0.5, minor=True)
            ax_heat.grid(which="minor", color="white", linewidth=0.8)
            ax_heat.tick_params(which="minor", bottom=False, left=False)
            for i in range(len(ids)):
                for j in range(len(conditions)):
                    v = matrix[i, j]
                    label = "NA" if not math.isfinite(v) else ("<0.1" if 0 < v < 0.1 else f"{v:.1f}")
                    ax_heat.text(j, i, label,
                                 ha="center", va="center", fontsize=7,
                                 color="white" if math.isfinite(v) and v > 57 else "#17232e")
            for ax in (ax_map, ax_heat):
                for spine in ax.spines.values():
                    spine.set_visible(False)
            bar_ax = fig.add_axes([0.75, 0.12, 0.19, 0.015])
            fig.colorbar(im, cax=bar_ax, orientation="horizontal", ticks=[0, 50, 100])
            bar_ax.tick_params(labelsize=7, length=2)
            bar_ax.set_xlabel("Reads within this deletion length (%)", fontsize=7)
            cov = {r["condition"]: r for r in coverage if int(r["length"]) == length and r["time"] == time}
            heat_position = ax_heat.get_position()
            for j, condition in enumerate(conditions):
                v = 100 * number(cov.get(condition, {}).get("top5_mean_fraction"))
                xpos = heat_position.x0 + heat_position.width * (j + 0.5) / len(conditions)
                fig.text(xpos, 0.19, "NA" if not math.isfinite(v) else f"{v:.1f}",
                         ha="center", va="center", fontsize=7)
            fig.text(0.11, 0.19, "Condition-specific Top5 coverage (%)",
                     ha="left", va="center", fontsize=7)
            fig.suptitle(f"{length} bp deletion patterns — {time.replace('h', ' h')}",
                         x=0.11, y=0.98, ha="left", fontsize=12)
            caption = (f"Cells: mean of independent sample fractions; denominator: all pure {length} bp DELETION reads.\n"
                       "Rows: union of condition-specific Top5; unselected cells retain their actual frequencies.\n"
                       "† Internal-only deletion; MH = NA. Coordinates are equivalent representatives, not unique cut sites.")
            fig.text(0.11, 0.015, caption, fontsize=7, va="bottom", linespacing=1.4)
            stem = f"modify_top5_{length}bp_{time}"
            save(fig, stem, args)
            (args.outdir / f"{stem}_caption.txt").write_text(caption + "\n", encoding="utf-8")
            with (args.outdir / f"{stem}_caption.txt").open("a", encoding="utf-8") as handle:
                handle.write("MH is annotated on the SIQ-reconstructed allele. Raw-mismatch reads remain in these pattern counts. "
                             "Displayed percentages are rounded; positive values below 0.1% are labelled <0.1.\n")
            created.append(stem)
    return created


def plot_mh(summary_dir, args):
    rows = read_rows(summary_dir / "exact_length_mh_by_group.csv")
    conditions = conditions_for(rows, args.conditions)
    names = label_map(rows)
    metric = f"mh_ge{args.mh_threshold}_fraction"
    if rows and metric + "_mean" not in rows[0]:
        raise ValueError(f"Missing expected MH column {metric}_mean")
    fig, axes = plt.subplots(len(args.lengths), len(args.times), squeeze=False,
                             figsize=(args.width_mm / 25.4, 2.45 * len(args.lengths) + 0.6),
                             sharey=True)
    lookup = {(int(r["deletion_length_bp"]), r["time"], r["condition"]): r for r in rows}
    colors = ["#b84779" if c == "Cy5" else "#507f9d" for c in conditions]
    limits = [0, 100]
    for key, row in lookup.items():
        if key[0] in args.lengths and key[1] in args.times and key[2] in conditions:
            m = 100 * number(row.get(metric + "_mean"))
            sd = 100 * number(row.get(metric + "_sd"))
            if math.isfinite(m):
                spread = sd if math.isfinite(sd) else 0
                limits.extend([m - spread, m + spread])
    for i, length in enumerate(args.lengths):
        for j, time in enumerate(args.times):
            ax = axes[i, j]
            selected = [lookup.get((length, time, c), {}) for c in conditions]
            means = np.array([100 * number(r.get(metric + "_mean")) for r in selected])
            sds = np.array([100 * number(r.get(metric + "_sd")) for r in selected])
            x = np.arange(len(conditions))
            ax.bar(x, means, color=colors, width=0.70, zorder=2)
            valid = np.isfinite(means) & np.isfinite(sds)
            ax.errorbar(x[valid], means[valid], yerr=sds[valid], fmt="none", color="#222222",
                        capsize=2, linewidth=0.8, zorder=3)
            for k, value in enumerate(means):
                if not math.isfinite(value):
                    ax.text(k, 2, "NA", ha="center", fontsize=7)
            ax.set_xticks(x, [names[c] for c in conditions], rotation=55, ha="right", fontsize=8)
            ax.set_title(f"{length} bp deletion, {time.replace('h', ' h')}", fontsize=10, loc="left")
            ax.set_ylim(min(0, min(limits) - 3), max(105, max(limits) + 3))
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.tick_params(axis="y", labelsize=8)
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color="#e4e7e9", linewidth=0.6)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if j == 0:
                ax.set_ylabel(f"MH candidate ≥ {args.mh_threshold} bp (%)", fontsize=9)
    fig.suptitle("Junction-compatible microhomology at matched deletion lengths", fontsize=12, y=0.99)
    caption = ("Mean ± sample SD across independent biological replicates; samples contribute equal weight.\n"
               "Denominator: strict MH-evaluable pure DELETION reads at each length (full Raw match, primers preserved, junction-compatible).\n"
               "Sequence-compatible MH candidates do not demonstrate use of a specific repair pathway.")
    fig.subplots_adjust(left=0.09, right=0.99, top=0.91, bottom=0.22, hspace=0.80, wspace=0.18)
    fig.text(0.09, 0.015, caption, fontsize=7, va="bottom", linespacing=1.4)
    stem = f"modify_exact_length_mh_ge{args.mh_threshold}"
    save(fig, stem, args)
    (args.outdir / f"{stem}_caption.txt").write_text(caption + "\n", encoding="utf-8")
    return [stem]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Root output folder containing modify/summary and modify/patterns")
    parser.add_argument("--outdir", type=Path, required=True, help="New or empty figure directory")
    parser.add_argument("--lengths", type=int, nargs="+", default=[6, 12])
    parser.add_argument("--times", nargs="+", default=["12h", "24h"])
    parser.add_argument("--conditions", nargs="+", help="Optional condition aliases in display order")
    parser.add_argument("--mh-threshold", type=int, choices=[1, 2, 3], default=2)
    parser.add_argument("--width-mm", type=float, default=220)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--formats", nargs="+", choices=["svg", "pdf", "png"], default=["svg", "pdf", "png"])
    parser.add_argument("--cmap", default="Blues")
    args = parser.parse_args(argv)
    if args.outdir.exists() and any(args.outdir.iterdir()):
        parser.error("Output directory is not empty; choose a new --outdir.")
    if args.width_mm < 140 or args.dpi < 72 or any(x <= 0 for x in args.lengths):
        parser.error("Use width >= 140 mm, dpi >= 72 and positive deletion lengths.")
    args.outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.unicode_minus": True, "savefig.transparent": False})
    created = []
    summary = args.results / "modify" / "summary"
    patterns = args.results / "modify" / "patterns"
    if (summary / "exact_length_mh_by_group.csv").is_file():
        created += plot_mh(summary, args)
    if (patterns / "length_group_patterns.csv").is_file():
        created += plot_top5(patterns, args)
    if not created:
        parser.error("No supported modification summary/pattern tables were found.")
    metadata = {"figure_stems": created, "formats": args.formats, "dpi": args.dpi,
                "width_mm_before_tight_crop": args.width_mm, "lengths": args.lengths,
                "times": args.times, "mh_threshold": args.mh_threshold,
                "matplotlib_version": matplotlib.__version__, "numpy_version": np.__version__,
                "scope": "Static numerical reproduction; not the Figure Studio UI or pixel-identical historical layouts.",
                "heatmap_scale_percent": [0, 100]}
    (args.outdir / "figure_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Created {len(created)} figures in {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
