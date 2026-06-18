#!/usr/bin/env python3
"""Materialize Orb+PLUMED metadynamics production runs per run dir.

Translation of qmmd's meta.py (meta-prep / meta-submit) to the ORB/ASE stack,
sharing the nested YAML vocabulary of plumed_cv.py
(method_dir / cv_dir / traj_name / run_ids).

Per run-<i> dir it:
  - requires <cv_dir>/plumed.dat (authored by `mlip plumed`);
  - asserts md.dt_fs == cv.dt_fs from cv_spec.yaml (PACE is baked into plumed.dat);
  - extracts the LAST frame of <replica_dir>/<traj_name> as start.xyz
    (extxyz carries Lattice/pbc, so the cell round-trips);
  - writes a production main.py wrapping ORBCalculator in ase ... Plumed;
  - reuses mlip.orb's run.sh / slurm.sh writers.
"""

from __future__ import annotations

import sys
import yaml
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

from ase.io import read, write

from mlip.orb import (
    RuntimeConfig,
    SlurmConfig,
    SlurmJobConfig,
    find_repo_root,
    orb_root_dir,
    write_run_sh,
    write_slurm_sh,
    submit_slurm,
)


@dataclass(frozen=True, slots=True)
class MetaPrepConfig:
    system: str
    buffer: float
    prefix: str
    method_dir: str
    bench_tag: str
    cv_dir: str
    replica_dir: str
    traj_name: str
    start_name: str
    run_ids: List[int]
    runtime: RuntimeConfig
    slurm: Optional[SlurmConfig] = None


def parse_range(value: Union[str, List[int], Tuple[int, int]]) -> Tuple[int, int]:
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return int(value[0]), int(value[1])
    if isinstance(value, str) and "-" in value:
        a, b = value.split("-", 1)
        return int(a.strip()), int(b.strip())
    raise RuntimeError(f"Bad range spec: {value!r}")


def parse_run_ids(value: Any) -> List[int]:
    if isinstance(value, int):
        return [int(value)]
    if isinstance(value, str):
        if "-" in value:
            a, b = parse_range(value)
            return list(range(a, b + 1))
        return [int(value)]
    if isinstance(value, list):
        out: List[int] = []
        for item in value:
            out.extend(parse_run_ids(item))
        return out
    raise RuntimeError(f"Bad run_ids spec: {value!r}")


def load_config(yaml_path: Path) -> MetaPrepConfig:
    data = yaml.safe_load(yaml_path.read_text())
    slurm_cfg = None
    if data.get("slurm") is not None:
        slurm_cfg = SlurmConfig(job=SlurmJobConfig(**data["slurm"]["job"]))
    return MetaPrepConfig(
        system=data["system"],
        buffer=float(data["buffer"]),
        prefix=data.get("prefix", "solv"),
        method_dir=data.get("method_dir", "orb"),
        bench_tag=data["bench_tag"],
        cv_dir=data["cv_dir"],
        replica_dir=data.get("replica_dir", "equil"),
        traj_name=data.get("traj_name", "traj.xyz"),
        start_name=data.get("start_name", "start.xyz"),
        run_ids=parse_run_ids(data["run_ids"]),
        runtime=RuntimeConfig(**data["runtime"]),
        slurm=slurm_cfg,
    )


class _orb_view:
    """Shim so we can reuse mlip.orb.orb_root_dir (expects .orb_dirname)."""
    def __init__(self, cfg: MetaPrepConfig):
        self.system = cfg.system
        self.buffer = cfg.buffer
        self.prefix = cfg.prefix
        self.orb_dirname = cfg.method_dir


def bench_dir(cfg: MetaPrepConfig, repo_root: Path) -> Path:
    return orb_root_dir(_orb_view(cfg), repo_root) / cfg.bench_tag


def extract_start_frame(equil_dir: Path, traj_name: str, dst: Path) -> None:
    traj = equil_dir / traj_name
    if not traj.exists() or traj.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty equil trajectory: {traj}")
    atoms = read(traj, index=-1)
    write(dst, atoms, format="extxyz")  # preserves Lattice + pbc


def assert_dt_consistent(cv_dir: Path, prod_spec: dict) -> None:
    """PACE/HEIGHT in plumed.dat were built for cv.dt_fs; the production
    md.dt_fs must match or the deposition cadence is silently wrong."""
    cv_spec_path = cv_dir / "cv_spec.yaml"
    if not cv_spec_path.exists():
        return
    cv_data = yaml.safe_load(cv_spec_path.read_text()) or {}
    cv_dt = (cv_data.get("cv") or {}).get("dt_fs")          # nested under cv:
    md_dt = (prod_spec.get("md") or {}).get("dt_fs")
    if cv_dt is not None and md_dt is not None and float(cv_dt) != float(md_dt):
        raise RuntimeError(
            f"dt_fs mismatch in {cv_dir}: plumed.dat built for cv.dt_fs={cv_dt} "
            f"but production md.dt_fs={md_dt}. Re-run `mlip cv` with matching dt_fs."
        )


def write_meta_main_py(out_dir: Path) -> Path:
    py = out_dir / "main.py"
    py.write_text(_META_MAIN_TEMPLATE)
    py.chmod(0o755)
    return py


def run_orb_meta_prep(yaml_path: Path) -> None:
    cfg = load_config(yaml_path)
    repo_root = find_repo_root(yaml_path)
    bench = bench_dir(cfg, repo_root)

    yaml_text = yaml_path.read_text()
    prod_spec = yaml.safe_load(yaml_text)

    for run_id in cfg.run_ids:
        run_root = bench / f"run-{run_id}"
        equil_dir = run_root / cfg.replica_dir
        cv_dir = run_root / cfg.cv_dir

        if not (cv_dir / "plumed.dat").exists():
            raise RuntimeError(f"Missing plumed.dat in {cv_dir} (run `mlip cv` first)")
        if (cv_dir / "main.py").exists():
            print(f"SKIP: {cv_dir} already prepared; not touching")
            continue

        assert_dt_consistent(cv_dir, prod_spec)
        extract_start_frame(equil_dir, cfg.traj_name, cv_dir / cfg.start_name)

        (cv_dir / "spec.yaml").write_text(yaml_text)       # consumed by main.py
        (cv_dir / "meta_spec.yaml").write_text(yaml_text)  # submit-match guard
        write_meta_main_py(cv_dir)
        write_run_sh(cfg, cv_dir)
        write_slurm_sh(cfg, cv_dir)
        print(f"OK: prepared meta run in {cv_dir} (run {run_id})")


def run_orb_meta_submit(yaml_path: Path) -> None:
    cfg = load_config(yaml_path)
    repo_root = find_repo_root(yaml_path)
    bench = bench_dir(cfg, repo_root)
    yaml_text = yaml_path.read_text()

    if cfg.slurm is None:
        print("NOTE: slurm config not provided; no submissions made")
        return

    targets: List[Path] = []
    for run_id in cfg.run_ids:
        cv_dir = bench / f"run-{run_id}" / cfg.cv_dir
        spec = cv_dir / "meta_spec.yaml"
        if not spec.exists() or spec.read_text() != yaml_text:
            print(f"SKIP: {cv_dir} meta_spec.yaml mismatch (not submitting)")
            continue
        if not (cv_dir / "slurm.sh").exists():
            print(f"SKIP: missing slurm.sh in {cv_dir}")
            continue
        targets.append(cv_dir)

    if not targets:
        print("NOTE: no matching meta dirs found; nothing submitted")
        return

    print("Will submit the following meta dirs:")
    for t in targets:
        print(f"  - {t}")
    resp = input(f"Proceed to submit {len(targets)} jobs? [y/N] ").strip().lower()
    if resp not in ("y", "yes"):
        print("Cancelled by user.")
        return
    for cv_dir in targets:
        print(f"Submitting via sbatch for {cv_dir}...")
        submit_slurm(cv_dir / "slurm.sh")
        print("OK: job submitted")


# --------------------------------------------------------------------------- #
# production driver template (written into each cv_dir as main.py)
# --------------------------------------------------------------------------- #
_META_MAIN_TEMPLATE = r'''#!/usr/bin/env python3
import secrets
from pathlib import Path

import numpy as np
import yaml

from ase import units
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.plumed import Plumed

import torch
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator


def setup_device(d):
    if d == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(d)


def main():
    cfg = yaml.safe_load(Path("spec.yaml").read_text())
    md, out, orb = cfg["md"], cfg["output"], cfg.get("orb", {})

    dt_fs = float(md["dt_fs"])
    T = float(md["T"])
    dt = dt_fs * units.fs
    friction = float(md["friction_per_fs"]) / units.fs
    n_prod = int(round(float(md["prod_ps"]) * 1000.0 / dt_fs))

    atoms = read("start.xyz", index=-1)  # extxyz carries Lattice + pbc
    if atoms.get_cell().rank == 3:
        atoms.set_pbc(True)
        atoms.wrap()
    atoms.info["charge"] = int(orb.get("charge", 0))
    atoms.info["spin"] = int(orb.get("spin", 1))

    model_name = orb.get("model", "orb_v3_conservative_omol")
    device = setup_device(orb.get("device", "auto"))
    orbff = getattr(pretrained, model_name)(
        device=device,
        precision=orb.get("precision", "float32-high"),
        compile=bool(orb.get("compile", False)),
    )
    orb_calc = ORBCalculator(orbff, device=device)

    plumed_input = Path("plumed.dat").read_text().splitlines()
    atoms.calc = Plumed(
        calc=orb_calc,
        input=plumed_input,
        timestep=dt,
        atoms=atoms,
        kT=units.kB * T,
    )

    seed = secrets.randbelow(2**31 - 1) + 1
    print(f"Using RNG seed (velocities only): {seed}")
    rng = np.random.get_state()
    np.random.seed(seed)
    MaxwellBoltzmannDistribution(atoms, temperature_K=T)
    np.random.set_state(rng)

    dyn = Langevin(atoms, dt, temperature_K=T, friction=friction)
    log_stride = int(out["log_stride"])
    traj_stride = int(out["traj_stride"])

    with open("md.log", "w") as flog:
        flog.write("# step time_ps Epot_eV Ekin_eV T_K\n")

        def log_cb():
            if dyn.nsteps % log_stride:
                return
            tps = dyn.nsteps * dt_fs * 1e-3
            flog.write(
                f"{dyn.nsteps} {tps:.6f} "
                f"{atoms.get_potential_energy():.8f} "
                f"{atoms.get_kinetic_energy():.8f} "
                f"{atoms.get_temperature():.3f}\n"
            )
            flog.flush()

        def traj_cb():
            if dyn.nsteps % traj_stride == 0:
                if atoms.get_cell().volume > 0.0:
                    atoms.wrap()
                write("traj.xyz", atoms, append=True)

        dyn.attach(log_cb, interval=1)
        dyn.attach(traj_cb, interval=1)
        print(f"Production steps: {n_prod} (prod_ps={md['prod_ps']}, dt_fs={dt_fs})")
        dyn.run(n_prod)

    print("OK: wrote md.log, traj.xyz; PLUMED wrote HILLS, COLVAR")


if __name__ == "__main__":
    main()
'''


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python3 orb_meta.py configs/SYS/orb-meta/orb-meta.yaml", file=sys.stderr)
        raise SystemExit(2)
    run_orb_meta_prep(Path(sys.argv[1]).resolve())


if __name__ == "__main__":
    main()
