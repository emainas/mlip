#!/usr/bin/env python3
"""Compare Slurm wall time and TeraChem processing time across QM/MM frames."""

import argparse
import re
import statistics
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import yaml


def slurm_records(job_id: str) -> dict[int, dict[str, str]]:
    result = subprocess.run(
        [
            "sacct",
            "-j",
            job_id,
            "-X",
            "-n",
            "-P",
            "--format=JobID,ElapsedRaw,NodeList,AllocTRES,State",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    records: dict[int, dict[str, str]] = {}
    for line in result.stdout.splitlines():
        fields = line.split("|")
        match = re.fullmatch(rf"{re.escape(job_id)}_(\d+)", fields[0])
        if not match:
            continue
        task = int(match.group(1))
        records[task] = {
            "wall": fields[1],
            "node": fields[2],
            "resources": fields[3],
            "state": fields[4],
        }
    return records


def terachem_record(frame_dir: Path) -> dict[str, float | int | str]:
    text = (frame_dir / "tc.out").read_text()
    processing = re.search(r"Total processing time:\s+([0-9.]+) sec", text)
    device = re.search(r"Device 0:\s+([^,]+),\s+([0-9]+)MB,\s+CC\s+([0-9.]+)", text)
    aos = re.search(r"Total AOs:\s+(\d+)", text)
    iterations = [int(value) for value in re.findall(r"^\|\s+(\d+)\s+", text, re.MULTILINE)]
    if not processing or not device or not aos or not iterations:
        raise RuntimeError(f"Incomplete TeraChem timing metadata in {frame_dir}")
    manifest = yaml.safe_load((frame_dir / "manifest.yaml").read_text())
    qm_atoms = int(manifest["qm_atoms"])
    return {
        "frame": int(manifest["source_frame"]),
        "qm_atoms": qm_atoms,
        "mm_atoms": int(manifest["carved"]["atoms"]) - qm_atoms,
        "processing": float(processing.group(1)),
        "gpu": device.group(1).strip(),
        "gpu_memory_mb": int(device.group(2)),
        "compute_capability": device.group(3),
        "aos": int(aos.group(1)),
        "iterations": max(iterations),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("job_id")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    frame_dirs = sorted(args.input_dir.glob("frame_*"))
    slurm = slurm_records(args.job_id)
    if len(frame_dirs) != 10 or len(slurm) != 10:
        raise RuntimeError(f"Expected 10 frames and Slurm tasks; found {len(frame_dirs)} and {len(slurm)}")

    records = []
    for task, frame_dir in enumerate(frame_dirs):
        record = terachem_record(frame_dir)
        record.update(
            wall=float(slurm[task]["wall"]),
            node=slurm[task]["node"],
            resources=slurm[task]["resources"],
        )
        records.append(record)

    cmap = plt.get_cmap("viridis")
    colors = [cmap(index / (len(records) - 1)) for index in range(len(records))]
    positions = list(range(len(records)))
    wall = [float(record["wall"]) for record in records]
    processing = [float(record["processing"]) for record in records]

    fig, ax = plt.subplots(figsize=(16, 9))
    bars = ax.bar(positions, wall, width=0.72, color=colors, edgecolor="black", linewidth=1.0)
    ax.scatter(
        positions,
        processing,
        color="black",
        marker="_",
        s=650,
        linewidths=3.0,
        zorder=3,
        label="TeraChem processing time",
    )
    for bar, seconds in zip(bars, wall):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            seconds + 1.2,
            f"{seconds:.0f} s",
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    labels = [
        f"Frame {record['frame']}\n{record['qm_atoms']} QM atoms\n"
        f"{record['mm_atoms']:,} MM atoms\n{record['iterations']} SCF iter."
        for record in records
    ]
    ax.set_xticks(positions, labels)
    ax.set_ylabel("Wall-clock time (seconds)", fontsize=17)
    ax.set_xlabel("Trajectory frame", fontsize=17, labelpad=12)
    ax.tick_params(axis="x", labelsize=11)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylim(0, max(wall) * 1.28)
    ax.grid(False)

    first = records[0]
    resources = str(first["resources"])
    cpu = re.search(r"(?:^|,)cpu=(\d+)", resources).group(1)
    memory = re.search(r"(?:^|,)mem=([^,]+)", resources).group(1)
    node = str(first["node"])
    ax.set_title(
        f"{first['gpu']} ({first['gpu_memory_mb']} MB, CC {first['compute_capability']}) on {node}"
        f"  |  1 GPU, {cpu} CPUs, {memory} RAM per task\n"
        "B3LYP/def2-SVP  |  spherical basis  |  8 link atoms in addition to the reported QM-region atoms",
        fontsize=13,
        pad=16,
    )
    fig.suptitle("Osteocalcin QM/MM single-point wall-clock comparison", fontsize=21, y=0.98)

    overhead = [total - compute for total, compute in zip(wall, processing)]
    ax.text(
        0.5,
        0.96,
        f"Mean wall time: {statistics.mean(wall):.1f} s\n"
        f"Mean TeraChem time: {statistics.mean(processing):.1f} s\n"
        f"Mean launch/exit overhead: {statistics.mean(overhead):.1f} s",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=12,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "white", "edgecolor": "0.7", "alpha": 0.9},
    )
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="black",
                marker="_",
                linestyle="None",
                markersize=22,
                markeredgewidth=3,
                label="TeraChem processing time",
            )
        ],
        loc="upper right",
        fontsize=12,
        frameon=False,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
