#!/usr/bin/env python3
"""Plot dihedral time series with per-dihedral subplots and overlaid runs.

Applies a running average to each run (default window=50).
"""

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt


def read_series(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise RuntimeError(f"Bad dihedral file: {path}")
    return data[:, 1]


def _dih_filename(label: str) -> str:
    label = label.strip()
    if label.startswith("dih_") and label.endswith(".dat"):
        return label
    if label.startswith("dih_"):
        return f"{label}.dat"
    return f"dih_{label}.dat"


def discover_runs(runs_path: Path, rel_dir: Path, labels: List[str]) -> List[Tuple[int, List[Path]]]:
    out: List[Tuple[int, List[Path]]] = []
    for run_dir in sorted(runs_path.glob("run-*")):
        if not run_dir.is_dir():
            continue
        try:
            run_id = int(run_dir.name.split("-")[-1])
        except Exception:
            continue
        analysis_dir = run_dir / rel_dir
        paths: List[Path] = []
        missing = False
        for label in labels:
            p = analysis_dir / _dih_filename(label)
            if not p.exists():
                print(f"WARN: missing {p.name} in {analysis_dir}")
                missing = True
                break
            paths.append(p)
        if missing:
            continue
        out.append((run_id, paths))
    return out


def running_average(y: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return y
    w = int(window)
    if w % 2 == 0:
        w += 1
    kernel = np.ones(w) / w
    return np.convolve(y, kernel, mode="same")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-path", required=True, type=Path, help="Path containing run-* directories")
    ap.add_argument("--rel-dir", default="equil/analysis", help="Relative path under each run")
    ap.add_argument("--labels", default="chi1,chi2",
                    help="Comma-separated dihedral labels (e.g., chi1,chi2 or single5,single10,single15)")
    ap.add_argument("--out", type=Path, default=Path("reports/dihedral_timeseries.png"))
    ap.add_argument("--window", type=int, default=50, help="Running average window (frames)")
    ap.add_argument("--style", type=Path, default=Path("plotting/prl.mplstyle"),
                    help="Matplotlib style file")
    args = ap.parse_args()

    labels = [x.strip() for x in args.labels.split(",") if x.strip()]
    if not labels:
        raise SystemExit("No labels provided")

    if args.style.exists():
        plt.style.use(args.style)

    runs = discover_runs(args.runs_path, Path(args.rel_dir), labels)
    if not runs:
        raise SystemExit("No runs found")

    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharey=True, sharex=True)
    axes_list = axes.flatten()

    cmap = plt.get_cmap("tab20")
    for idx, (run_id, paths) in enumerate(runs):
        color = cmap(idx % cmap.N)
        for ax, label, path in zip(axes_list, labels, paths):
            series = read_series(path)
            series = running_average(series, args.window)
            t = np.arange(len(series))
            ax.scatter(t, series, s=6.0, color=color, alpha=0.7, label=f"run-{run_id}")
            ax.set_title(label)
            ax.grid(True, alpha=0.3)

    for i in range(len(labels), len(axes_list)):
        axes_list[i].set_axis_off()

    axes_list[0].set_ylabel("dihedral (deg)")
    axes_list[3].set_ylabel("dihedral (deg)")
    axes_list[-1].set_xlabel("frame")

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)


if __name__ == "__main__":
    main()
