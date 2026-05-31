#!/usr/bin/env python3
"""Plot Euler3D profiling CSV files with matplotlib."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TOP_LEVEL = [
    ("Advect velocity", "advect_velocity_ms"),
    ("Velocity source", "velocity_source_ms"),
    ("Boundary", "boundary_conditions_ms"),
    ("Projection", "projection_ms"),
    ("Advect scalars", "advect_scalars_ms"),
    ("Scalar source", "scalar_source_ms"),
]

PROJECTION = [
    ("Divergence", "projection_divergence_ms"),
    ("Build system", "pressure_system_build_ms"),
    ("Pressure solve", "pressure_solve_ms"),
    ("Apply gradient", "pressure_apply_gradient_ms"),
    ("Boundary", "projection_boundary_ms"),
    ("Validation", "projection_validation_ms"),
]

SCALAR = [
    ("Scalar advection", "scalar_advection_ms"),
    ("Postprocess", "scalar_postprocess_ms"),
]

TIME_SERIES = [
    ("Total", "total_ms"),
    ("Advect velocity", "advect_velocity_ms"),
    ("Projection", "projection_ms"),
    ("Advect scalars", "advect_scalars_ms"),
    ("Pressure solve", "pressure_solve_ms"),
]

STACKED = [
    ("Advect velocity", "advect_velocity_ms"),
    ("Velocity source", "velocity_source_ms"),
    ("Boundary", "boundary_conditions_ms"),
    ("Projection", "projection_ms"),
    ("Advect scalars", "advect_scalars_ms"),
    ("Scalar source", "scalar_source_ms"),
]


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items() if value != ""}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise ValueError(f"No profiling rows found in {path}")
    return rows


def read_summary(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows: list[dict[str, float | str]] = []
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "stage": row["stage"],
                    "total_ms": float(row["total_ms"]),
                    "ms_per_frame": float(row["ms_per_frame"]),
                    "percent_total": float(row["percent_total"]),
                }
            )
    if not rows:
        raise ValueError(f"No summary rows found in {path}")
    return rows


def average(rows: list[dict[str, float]], key: str) -> float:
    return sum(row.get(key, 0.0) for row in rows) / len(rows)


def save_average_bar(
    rows: list[dict[str, float]],
    stages: list[tuple[str, str]],
    baseline_key: str,
    title: str,
    output_path: Path,
) -> None:
    labels = [name for name, _ in stages]
    values = [average(rows, key) for _, key in stages]
    baseline = max(average(rows, baseline_key), 1e-12)
    colors = plt.cm.Set2(range(len(labels)))

    fig, ax = plt.subplots(figsize=(11, 5.8), layout="constrained")
    bars = ax.barh(labels, values, color=colors)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=15, weight="bold")
    ax.set_xlabel("Average time per frame (ms)")
    ax.grid(axis="x", alpha=0.25)

    for bar, value in zip(bars, values):
        percent = 100.0 * value / baseline
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2.0,
            f"  {value:.3f} ms ({percent:.1f}%)",
            va="center",
            fontsize=9,
        )

    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def save_time_share(rows: list[dict[str, float]], output_path: Path) -> None:
    labels = [name for name, _ in TOP_LEVEL]
    values = [average(rows, key) for _, key in TOP_LEVEL]
    total = max(sum(values), 1e-12)
    pie_labels = [label if value / total >= 0.01 else "" for label, value in zip(labels, values)]
    fig, (pie_ax, bar_ax) = plt.subplots(1, 2, figsize=(13, 5.8), layout="constrained")

    pie_ax.pie(
        values,
        labels=pie_labels,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 1.0 else "",
        startangle=90,
        counterclock=False,
    )
    pie_ax.set_title("Average time share", fontsize=14, weight="bold")

    bars = bar_ax.barh(labels, values, color=plt.cm.Set2(range(len(labels))))
    bar_ax.invert_yaxis()
    bar_ax.set_title("Average ms/frame", fontsize=14, weight="bold")
    bar_ax.set_xlabel("ms")
    bar_ax.grid(axis="x", alpha=0.25)
    for bar, value in zip(bars, values):
        bar_ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2.0, f" {value:.3f}", va="center")

    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def save_time_series(rows: list[dict[str, float]], output_path: Path) -> None:
    frames = [int(row["frame"]) for row in rows]
    fig, ax = plt.subplots(figsize=(12, 6), layout="constrained")
    for label, key in TIME_SERIES:
        ax.plot(frames, [row.get(key, 0.0) for row in rows], marker="o", linewidth=2, label=label)
    ax.set_title("Stage time over measured frames", fontsize=15, weight="bold")
    ax.set_xlabel("Measured frame")
    ax.set_ylabel("Time (ms)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def save_stacked_area(rows: list[dict[str, float]], output_path: Path) -> None:
    frames = [int(row["frame"]) for row in rows]
    labels = [label for label, _ in STACKED]
    values = [[row.get(key, 0.0) for row in rows] for _, key in STACKED]

    fig, ax = plt.subplots(figsize=(12, 6), layout="constrained")
    ax.stackplot(frames, values, labels=labels, alpha=0.9)
    ax.set_title("Stacked stage time per frame", fontsize=15, weight="bold")
    ax.set_xlabel("Measured frame")
    ax.set_ylabel("Time (ms)")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper left")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def save_summary(summary_rows: list[dict[str, float | str]], output_path: Path) -> None:
    rows = [row for row in summary_rows if row["stage"] != "total step"]
    rows = sorted(rows, key=lambda row: float(row["ms_per_frame"]))
    labels = [str(row["stage"]) for row in rows]
    values = [float(row["ms_per_frame"]) for row in rows]
    percents = [float(row["percent_total"]) for row in rows]

    fig_height = max(6.0, 0.36 * len(rows) + 1.8)
    fig, ax = plt.subplots(figsize=(12, fig_height), layout="constrained")
    bars = ax.barh(labels, values, color=plt.cm.tab20(range(len(rows))))
    ax.set_title("Euler3D profiling summary", fontsize=15, weight="bold")
    ax.set_xlabel("Average time per frame (ms)")
    ax.grid(axis="x", alpha=0.25)

    for bar, value, percent in zip(bars, values, percents):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2.0,
            f"  {value:.3f} ms ({percent:.1f}%)",
            va="center",
            fontsize=9,
        )

    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Euler3D profiling CSV files.")
    parser.add_argument("--input", default="outputs/euler3d/euler3d_profile_frames.csv", type=Path)
    parser.add_argument("--summary", default=None, type=Path)
    parser.add_argument("--out", default="outputs/euler3d", type=Path)
    args = parser.parse_args()

    rows = read_rows(args.input)
    summary_path = args.summary or args.input.with_name("euler3d_profile_summary.csv")
    summary_rows = read_summary(summary_path) if summary_path.exists() else None
    args.out.mkdir(parents=True, exist_ok=True)

    outputs = [
        args.out / "euler3d_time_share.png",
        args.out / "euler3d_projection_breakdown.png",
        args.out / "euler3d_scalar_breakdown.png",
        args.out / "euler3d_time_series.png",
        args.out / "euler3d_stacked_stage_times.png",
    ]
    save_time_share(rows, outputs[0])
    save_average_bar(rows, PROJECTION, "projection_ms", "Projection breakdown", outputs[1])
    save_average_bar(rows, SCALAR, "advect_scalars_ms", "Scalar advection breakdown", outputs[2])
    save_time_series(rows, outputs[3])
    save_stacked_area(rows, outputs[4])
    if summary_rows is not None:
        outputs.append(args.out / "euler3d_summary.png")
        save_summary(summary_rows, outputs[-1])

    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
