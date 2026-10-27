#!/usr/bin/env python3
"""Prepare non-periodic TeraChem/OpenMM QM/MM snapshots from Amber MD."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import yaml


def find_repo_root(start: Path) -> Path:
    path = start.resolve()
    for parent in [path] + list(path.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    raise RuntimeError(f"Could not find repository root from {start}")


def _resolve(repo_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (repo_root / path).resolve()


def _require_tool(name: str) -> Path:
    executable = shutil.which(name)
    if not executable:
        raise RuntimeError(f"{name} was not found in PATH; load AmberTools first")
    return Path(executable).resolve()


def _amber_python() -> Path:
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        raise RuntimeError("AMBERHOME is not set; load AmberTools first")
    python = Path(amberhome) / "miniconda" / "bin" / "python"
    if not python.exists():
        raise RuntimeError(f"Amber Python was not found: {python}")
    return python


def _run(command: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        rendered = " ".join(command)
        raise RuntimeError(f"Command failed: {rendered}\n{exc.stdout or ''}") from exc


def trajectory_frame_count(cpptraj: Path, topology: Path, trajectory: Path) -> int:
    result = _run([str(cpptraj), "-p", str(topology), "-y", str(trajectory), "-tl"])
    match = re.search(r"Frames:\s*(\d+)", result.stdout)
    if not match:
        raise RuntimeError(f"Could not determine trajectory length:\n{result.stdout}")
    return int(match.group(1))


def evenly_spaced_frames(start: int, stop: int, count: int) -> list[int]:
    if start < 1 or stop < start:
        raise ValueError("Frame selection requires 1 <= start <= stop")
    if count < 1:
        raise ValueError("Frame selection count must be positive")
    if count == 1:
        return [start]
    frames = [round(start + index * (stop - start) / (count - 1)) for index in range(count)]
    if len(set(frames)) != count:
        raise ValueError("Requested more evenly spaced frames than the selected interval contains")
    return frames


def _residue_range(value: Any, key: str) -> tuple[int, int]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{key} must be a two-integer list")
    first, last = int(value[0]), int(value[1])
    if first < 1 or last < first:
        raise ValueError(f"Invalid {key}: {value}")
    return first, last


def _parse_pdb(path: Path) -> dict[int, list[tuple[str, str, float, float, float]]]:
    residues: dict[int, list[tuple[str, str, float, float, float]]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atom = line[12:16].strip()
            residue_name = line[17:20].strip()
            residue_number = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError as exc:
            raise RuntimeError(f"Could not parse generated PDB line: {line}") from exc
        residues.setdefault(residue_number, []).append((atom, residue_name, *xyz))
    if not residues:
        raise RuntimeError(f"No atoms were read from {path}")
    return residues


def _center_of_residues(
    residues: dict[int, list[tuple[str, str, float, float, float]]], numbers: list[int]
) -> tuple[float, float, float]:
    points = [atom[2:] for number in numbers for atom in residues.get(number, [])]
    if len(points) != len(numbers):
        raise RuntimeError(
            f"Expected one center atom in each of residues {numbers}; found {len(points)} atoms"
        )
    return tuple(sum(point[axis] for point in points) / len(points) for axis in range(3))


def _distance(point: tuple[float, float, float], center: tuple[float, float, float]) -> float:
    return math.sqrt(sum((point[index] - center[index]) ** 2 for index in range(3)))


def _coordinating_waters(
    residues: dict[int, list[tuple[str, str, float, float, float]]],
    calcium_residues: list[int],
    config: dict[str, Any],
) -> list[int]:
    cutoff = float(config["cutoff_angstrom"])
    if cutoff <= 0:
        raise ValueError("qm_region.coordinating_waters.cutoff_angstrom must be positive")
    calcium_points = []
    for number in calcium_residues:
        atoms = residues.get(number, [])
        points = [atom[2:] for atom in atoms if atom[0] == "CA" and atom[1] == "CA"]
        if len(points) != 1:
            raise RuntimeError(f"Expected one Ca atom in residue {number}; found {len(points)}")
        calcium_points.append(points[0])

    residue_name = str(config.get("residue_name", "WAT"))
    oxygen_name = str(config.get("oxygen_atom", "O"))
    selected = []
    for number, atoms in residues.items():
        oxygens = [
            atom[2:]
            for atom in atoms
            if atom[0] == oxygen_name and atom[1] == residue_name
        ]
        if any(_distance(oxygen, calcium) <= cutoff for oxygen in oxygens for calcium in calcium_points):
            selected.append(number)
    return sorted(selected)


def _qm_mask(region: dict[str, Any]) -> str:
    calcium = ",".join(str(value) for value in region["calcium_residues"])
    selections = [f":{calcium}"]
    for residue_name, numbers in region["carboxylate_residues"].items():
        atoms = region["sidechain_atoms"].get(residue_name)
        if not atoms:
            raise ValueError(f"No sidechain atom list provided for {residue_name}")
        residues = ",".join(str(value) for value in numbers)
        selections.append(f":{residues}@{','.join(atoms)}")
    return "|".join(f"({selection})" for selection in selections)


def _selected_atoms(cpptraj: Path, topology: Path, mask: str) -> list[int]:
    result = _run([str(cpptraj), "-p", str(topology), "-ms", mask])
    match = re.search(r"Selected=\s*([\d\s]+)", result.stdout)
    if not match:
        raise RuntimeError(f"Amber mask selected no atoms: {mask}\n{result.stdout}")
    return [int(value) for value in match.group(1).split()]


def _write_terachem_input(path: Path, config: dict[str, Any]) -> None:
    ordered = [
        ("prmtop", "system.parm7"),
        ("coordinates", "frame.rst7"),
        ("qmindices", "region.qm"),
        ("basis", config["basis"]),
        ("method", config["method"]),
        ("charge", config["charge"]),
        ("spinmult", config["spinmult"]),
    ]
    for key in (
        "precision",
        "dftgrid",
        "dynamicgrid",
        "sphericalbasis",
        "threall",
        "convthre",
        "maxit",
    ):
        if key in config:
            ordered.append((key, config[key]))
    ordered.append(("run", config.get("run", "gradient")))
    for key, value in config.get("extra", {}).items():
        ordered.append((str(key), value))
    path.write_text("".join(f"{key:<16} {value}\n" for key, value in ordered) + "end\n")


def _extract_frame(
    cpptraj: Path,
    ambpdb: Path,
    topology: Path,
    trajectory: Path,
    frame: int,
    center_mask: str,
    frame_dir: Path,
) -> Path:
    script = frame_dir / "extract.cpptraj"
    full_restart = frame_dir / "full-frame.rst7"
    script.write_text(
        f"parm {topology}\n"
        f"trajin {trajectory} {frame} {frame}\n"
        f"autoimage anchor {center_mask}\n"
        f"center {center_mask} origin mass\n"
        "image origin center familiar\n"
        f"trajout {full_restart.name} restart nobox\n"
        "run\nquit\n"
    )
    result = _run([str(cpptraj), "-i", script.name], cwd=frame_dir)
    (frame_dir / "extract.cpptraj.out").write_text(result.stdout)
    pdb_result = _run(
        [str(ambpdb), "-p", str(topology), "-c", full_restart.name], cwd=frame_dir
    )
    pdb = frame_dir / "full-frame.pdb"
    pdb.write_text(pdb_result.stdout)
    return full_restart


def _worker(command: list[str]) -> dict[str, Any]:
    worker = Path(__file__).with_name("qmmm_worker.py").resolve()
    amberhome = Path(os.environ["AMBERHOME"]).resolve()
    site_packages = sorted((amberhome / "lib").glob("python*/site-packages"))
    if not site_packages:
        raise RuntimeError(f"ParmEd site-packages were not found under {amberhome / 'lib'}")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(site_packages[-1])
    try:
        result = subprocess.run(
            [str(_amber_python()), str(worker), *command],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"ParmEd worker failed:\n{exc.stdout or ''}") from exc
    try:
        return json.loads(result.stdout.splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not parse ParmEd worker output:\n{result.stdout}") from exc


def _write_runner(output_dir: Path, environment: str, executable: str) -> None:
    runner = output_dir / "run_all.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        f"source \"{environment}\"\n\n"
        "set -euo pipefail\n\n"
        "while read -r frame; do\n"
        "  frame_dir=$(printf 'frame_%04d' \"$frame\")\n"
        "  echo \"==> $frame_dir\"\n"
        f"  (cd \"$frame_dir\" && \"{executable}\" tc.in > tc.out)\n"
        "done < frames.txt\n"
    )
    runner.chmod(0o755)


def _write_slurm_array(
    output_dir: Path,
    environment: str,
    executable: str,
    frame_count: int,
    config: dict[str, Any],
) -> None:
    max_parallel = int(config.get("max_parallel", frame_count))
    array = f"0-{frame_count - 1}%{max_parallel}"
    mem_directive = f"#SBATCH --mem={config['mem']}\n" if config.get("mem") else ""
    script = output_dir / "slurm_array.sh"
    script.write_text(
        "#!/usr/bin/env bash\n"
        f"#SBATCH -J {config.get('job_name', 'mlip-qmmm')}\n"
        f"#SBATCH -p {config['partition']}\n"
        "#SBATCH -N 1\n"
        "#SBATCH -n 1\n"
        f"#SBATCH -t {config['time']}\n"
        f"#SBATCH --qos={config['qos']}\n"
        f"#SBATCH --gres={config.get('gres', 'gpu:1')}\n"
        f"#SBATCH --cpus-per-task={int(config.get('cpus_per_task', 2))}\n"
        f"{mem_directive}"
        f"#SBATCH --array={array}\n"
        "#SBATCH -o slurm.%A_%a.out\n"
        "#SBATCH -e slurm.%A_%a.err\n\n"
        f"source \"{environment}\"\n"
        "set -euo pipefail\n"
        "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}\n"
        f"output_dir=\"{output_dir}\"\n"
        "mapfile -t frames < \"$output_dir/frames.txt\"\n"
        "frame=${frames[$SLURM_ARRAY_TASK_ID]}\n"
        "frame_dir=$(printf 'frame_%04d' \"$frame\")\n"
        "cd \"$output_dir/$frame_dir\"\n"
        f"\"{executable}\" tc.in > tc.out\n"
    )
    script.chmod(0o755)


def run_qmmm(yaml_path: Path) -> Path:
    yaml_path = yaml_path.resolve()
    data = yaml.safe_load(yaml_path.read_text())
    repo_root = find_repo_root(yaml_path)

    known_limitations = [str(note) for note in data.get("known_limitations", [])]
    for note in known_limitations:
        print(f"WARNING: {note}")

    mother = data["mother_run"]
    topology = _resolve(repo_root, mother["topology"])
    trajectory = _resolve(repo_root, mother["trajectory"])
    output_dir = _resolve(repo_root, data["output_dir"])
    for path in (topology, trajectory):
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")

    cpptraj = _require_tool("cpptraj")
    ambpdb = _require_tool("ambpdb")
    nframes = trajectory_frame_count(cpptraj, topology, trajectory)
    expected = int(mother.get("expected_frames", nframes))
    if nframes != expected:
        raise RuntimeError(
            f"Mother trajectory has {nframes} frames, expected {expected}. "
            "Wait for production to finish or correct mother_run.expected_frames."
        )

    selection = mother["selection"]
    if selection.get("mode") != "evenly_spaced":
        raise ValueError("Only mother_run.selection.mode: evenly_spaced is currently supported")
    frames = evenly_spaced_frames(
        int(selection.get("start", 1)),
        nframes if selection.get("stop", "last") == "last" else int(selection["stop"]),
        int(selection["count"]),
    )

    sphere = data["sphere"]
    radius = float(sphere["radius_angstrom"])
    if radius <= 0:
        raise ValueError("sphere.radius_angstrom must be positive")
    always_first, always_last = _residue_range(
        sphere["always_keep_residue_range"], "sphere.always_keep_residue_range"
    )
    center_residues = [int(value) for value in sphere["center_residues"]]
    center_mask = ":" + ",".join(str(value) for value in center_residues)

    qm_region = data["qm_region"]
    qm_mask = _qm_mask(qm_region)
    expected_base_qm_atoms = int(qm_region["expected_atoms"])
    expected_links = int(qm_region["expected_link_bonds"])
    water_config = qm_region.get("coordinating_waters")

    if output_dir.exists() and any(output_dir.iterdir()) and not data.get("overwrite", False):
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}. Set overwrite: true to replace frame files."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "frames.txt").write_text("".join(f"{frame}\n" for frame in frames))
    (output_dir / "spec.yaml").write_text(yaml_path.read_text())

    summaries = []
    for frame in frames:
        frame_dir = output_dir / f"frame_{frame:04d}"
        frame_dir.mkdir(parents=True, exist_ok=True)
        full_restart = _extract_frame(
            cpptraj, ambpdb, topology, trajectory, frame, center_mask, frame_dir
        )
        residues = _parse_pdb(frame_dir / "full-frame.pdb")
        center = _center_of_residues(residues, center_residues)
        qm_water_residues = (
            _coordinating_waters(residues, center_residues, water_config)
            if water_config
            else []
        )
        if water_config:
            minimum, maximum = _residue_range(
                water_config["expected_count_range"],
                "qm_region.coordinating_waters.expected_count_range",
            )
            if not minimum <= len(qm_water_residues) <= maximum:
                raise RuntimeError(
                    f"Frame {frame}: found {len(qm_water_residues)} coordinating waters; "
                    f"expected {minimum}-{maximum}"
                )

        kept = set(range(always_first, always_last + 1))
        for residue_number, atoms in residues.items():
            if always_first <= residue_number <= always_last:
                continue
            if min(_distance(atom[2:], center) for atom in atoms) <= radius:
                kept.add(residue_number)
        keep_file = frame_dir / "kept-residues.txt"
        keep_file.write_text("".join(f"{number}\n" for number in sorted(kept)))
        qm_water_file = frame_dir / "qm-water-residues.txt"
        qm_water_file.write_text("".join(f"{number}\n" for number in qm_water_residues))

        crop_summary = _worker(
            [
                "crop",
                "--topology",
                str(topology),
                "--coordinates",
                str(full_restart),
                "--keep-residues",
                str(keep_file),
                "--output-topology",
                str(frame_dir / "system.parm7"),
                "--output-coordinates",
                str(frame_dir / "frame.rst7"),
                *(
                    ["--tip3p-qm-water-residues", str(qm_water_file)]
                    if qm_water_residues
                    else []
                ),
                *(
                    [
                        "--ion-1264-mask",
                        str(data["nonbonded_1264"]["mask"]),
                        "--ion-1264-watermodel",
                        str(data["nonbonded_1264"]["watermodel"]),
                    ]
                    if data.get("nonbonded_1264")
                    else []
                ),
            ]
        )
        topology_text = (frame_dir / "system.parm7").read_text()
        if "%FLAG BOX_DIMENSIONS" in topology_text:
            raise RuntimeError(f"Periodic box flag survived carving for frame {frame}")
        if data.get("nonbonded_1264") and "%FLAG LENNARD_JONES_CCOEF" not in topology_text:
            raise RuntimeError(f"12-6-4 C coefficients are missing for frame {frame}")

        selected_one_based = _selected_atoms(cpptraj, frame_dir / "system.parm7", qm_mask)
        selected_zero_based = {index - 1 for index in selected_one_based}
        selected_zero_based.update(crop_summary["qm_water_atoms_zero_based"])
        expected_qm_atoms = expected_base_qm_atoms + 3 * len(qm_water_residues)
        if len(selected_zero_based) != expected_qm_atoms:
            raise RuntimeError(
                f"Frame {frame}: QM selection contains {len(selected_zero_based)} atoms; "
                f"expected {expected_qm_atoms}. Mask: {qm_mask}"
            )
        region = frame_dir / "region.qm"
        region.write_text("".join(f"{index}\n" for index in sorted(selected_zero_based)))
        inspect_summary = _worker(
            [
                "inspect",
                "--topology",
                str(frame_dir / "system.parm7"),
                "--coordinates",
                str(frame_dir / "frame.rst7"),
                "--qm-indices",
                str(region),
            ]
        )
        if len(inspect_summary["link_bonds"]) != expected_links:
            raise RuntimeError(
                f"Frame {frame}: found {len(inspect_summary['link_bonds'])} QM/MM link bonds; "
                f"expected {expected_links}"
            )

        tc = dict(data["terachem"])
        seed_base = tc.pop("seed_base", None)
        if seed_base is not None:
            tc.setdefault("extra", {})["seed"] = int(seed_base) + frame
        tc["charge"] = int(qm_region["charge"])
        tc["spinmult"] = int(qm_region["multiplicity"])
        _write_terachem_input(frame_dir / "tc.in", tc)
        manifest = {
            "source_frame": frame,
            "source_topology": str(topology),
            "source_trajectory": str(trajectory),
            "sphere_center_angstrom": list(center),
            "sphere_radius_angstrom": radius,
            "kept_residues": len(kept),
            "carved": crop_summary,
            "qm_mask": qm_mask,
            "qm_water_residues_source": qm_water_residues,
            "qm_water_count": len(qm_water_residues),
            "qm_charge": int(qm_region["charge"]),
            "qm_multiplicity": int(qm_region["multiplicity"]),
            "qm_atoms": len(selected_zero_based),
            "link_bonds": inspect_summary["link_bonds"],
        }
        (frame_dir / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
        summaries.append(manifest)
        print(
            f"frame {frame}: {crop_summary['atoms']} atoms, "
            f"{crop_summary['residues']} residues, {len(selected_zero_based)} QM atoms, "
            f"{len(qm_water_residues)} QM waters"
        )

    execution = data.get("execution", {})
    if execution:
        _write_runner(
            output_dir,
            str(_resolve(repo_root, execution["environment"])),
            str(_resolve(repo_root, execution["executable"])),
        )
        if execution.get("slurm"):
            _write_slurm_array(
                output_dir,
                str(_resolve(repo_root, execution["environment"])),
                str(_resolve(repo_root, execution["executable"])),
                len(frames),
                execution["slurm"],
            )
    (output_dir / "summary.yaml").write_text(
        yaml.safe_dump(
            {
                "frames": frames,
                "known_limitations": known_limitations,
                "snapshots": summaries,
            },
            sort_keys=False,
        )
    )
    print(f"OK: prepared {len(frames)} TeraChem QM/MM snapshots in {output_dir}")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="mlip qmmm", description="Pluck Amber MD frames into TeraChem QM/MM droplets"
    )
    parser.add_argument("yaml", type=Path)
    args = parser.parse_args()
    run_qmmm(args.yaml)


if __name__ == "__main__":
    main()
