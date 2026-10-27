#!/usr/bin/env python3
"""Plot TeraChem SCF energy convergence for prepared QM/MM frames."""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MaxNLocator
import yaml


def read_scf(path: Path) -> tuple[list[int], list[float], float]:
    steps: list[int] = []
    energies: list[float] = []
    final_energy = None
    for line in path.read_text().splitlines():
        if re.match(r"^\|\s+\d+\s+", line):
            fields = line.removeprefix("|").split()
            steps.append(int(fields[0]))
            energies.append(float(fields[5]))
        elif line.startswith("FINAL ENERGY:"):
            final_energy = float(line.split()[2])
    if not steps or final_energy is None:
        raise RuntimeError(f"Incomplete SCF output: {path}")
    return steps, energies, final_energy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--metric",
        choices=("energy", "energy-change"),
        default="energy",
        help="Plot absolute SCF energy or the log-scale absolute change per iteration",
    )
    parser.add_argument(
        "--overlay",
        action="store_true",
        help="Overlay all frames on one axis (intended for --metric energy-change)",
    )
    args = parser.parse_args()

    frame_dirs = sorted(args.input_dir.glob("frame_*"))
    if len(frame_dirs) != 10:
        raise RuntimeError(f"Expected 10 frame directories, found {len(frame_dirs)}")

    if args.overlay:
        if args.metric != "energy-change":
            raise ValueError("--overlay requires --metric energy-change")
        records = []
        for frame_dir in frame_dirs:
            steps, energies, _ = read_scf(frame_dir / "tc.out")
            manifest = yaml.safe_load((frame_dir / "manifest.yaml").read_text())
            records.append(
                (
                    int(manifest["source_frame"]),
                    int(manifest["qm_atoms"]),
                    steps[1:],
                    [abs(current - previous) for previous, current in zip(energies, energies[1:])],
                )
            )

        fig, ax = plt.subplots(figsize=(14, 9))
        norm = Normalize(vmin=records[0][0], vmax=records[-1][0])
        cmap = plt.get_cmap("viridis")
        for frame, qm_atoms, steps, changes in records:
            ax.plot(
                steps,
                changes,
                color=cmap(norm(frame)),
                marker="o",
                markersize=6.0,
                linewidth=3.0,
                alpha=0.9,
                label=f"Frame {frame} ({qm_atoms} QM atoms)",
            )
        ax.set_yscale("log")
        ax.set_xlim(1.5, 18.5)
        ax.xaxis.set_major_locator(FixedLocator(list(range(2, 19, 2))))
        ax.set_xlabel("SCF iteration", fontsize=17)
        ax.set_ylabel(r"Absolute energy change $|E_n-E_{n-1}|$ ($E_h$)", fontsize=17)
        ax.tick_params(axis="both", which="both", labelsize=15)
        ax.grid(False)
        ax.set_title(
            "Osteocalcin QM/MM SCF energy-change convergence\n"
            "B3LYP/def2-SVP, spherical basis",
            fontsize=21,
            pad=14,
        )
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.13),
            ncol=2,
            fontsize=11,
            frameon=False,
            handlelength=3.0,
        )
        fig.subplots_adjust(left=0.12, right=0.97, top=0.87, bottom=0.25)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    fig, axes = plt.subplots(5, 2, figsize=(15, 22), sharex=True, sharey=False)
    for ax, frame_dir in zip(axes.flat, frame_dirs):
        steps, energies, final_energy = read_scf(frame_dir / "tc.out")
        manifest = yaml.safe_load((frame_dir / "manifest.yaml").read_text())
        frame = int(manifest["source_frame"])
        qm_atoms = int(manifest["qm_atoms"])

        if args.metric == "energy-change":
            plot_steps = steps[1:]
            plot_values = [abs(current - previous) for previous, current in zip(energies, energies[1:])]
        else:
            plot_steps = steps
            plot_values = energies

        ax.plot(
            plot_steps,
            plot_values,
            color="#146C94",
            marker="o",
            markersize=6.5,
            linewidth=3.0,
        )
        ax.xaxis.set_major_locator(FixedLocator([1, 4, 7, 10, 13, 16, 18]))
        ax.set_xlim(0.5, 18.5)
        if args.metric == "energy-change":
            ax.set_yscale("log")
            ax.set_ylabel(r"$|E_n-E_{n-1}|$ ($E_h$)", fontsize=11)
            ax.grid(True, which="both", alpha=0.25, linewidth=0.7)
        else:
            ax.axhline(final_energy, color="black", linestyle="--", linewidth=1.5, alpha=0.75)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            span = max(energies) - min(energies)
            margin = max(0.05 * span, 0.05)
            ax.set_ylim(min(energies) - margin, max(energies) + margin)
            ax.set_ylabel(r"SCF energy ($E_h$)", fontsize=11)
            ax.grid(True, alpha=0.25, linewidth=0.7)
        ax.set_title(
            f"Trajectory frame {frame} ({qm_atoms} QM atoms)\n"
            rf"$E_{{final}}$ = {final_energy:.8f} $E_h$",
            fontsize=12,
            pad=8,
        )
        ax.set_xlabel("SCF iteration", fontsize=11)
        ax.tick_params(axis="both", which="both", labelsize=10, labelleft=True, labelbottom=True)

    title = (
        "Osteocalcin QM/MM SCF energy-change convergence\n"
        if args.metric == "energy-change"
        else "Osteocalcin QM/MM single-point SCF energy convergence\n"
    )
    fig.suptitle(
        title + "B3LYP/def2-SVP, spherical basis",
        fontsize=17,
        y=0.998,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=2.0, w_pad=1.6)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
