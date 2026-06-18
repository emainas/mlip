#!/usr/bin/env python3
"""Author a PLUMED metadynamics input (plumed.dat) per run dir.
 
Translation of qmmd's ncoord.py to the ORB/ASE/PLUMED stack.
 
Layout written:
  systems/<system>/<prefix>_<buffer>/<method_dir>/<bench_tag>/run-<i>/<cv_dir>/plumed.dat
 
Atom indices are 1-based (PLUMED convention) and pass through from the YAML
unchanged. Symbols are read from the last frame of <run-i>/<replica_dir>/<traj_name>
and the resolved element of every selected index is printed for confirmation.
"""
 
from __future__ import annotations
 
import sys
import yaml
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
 
from ase.io import read
 
from mlip.orb import find_repo_root, orb_root_dir
 
HARTREE_TO_EV = 27.211386245988
 
 
# --------------------------------------------------------------------------- #
# group selection (1-based output, PLUMED-ready)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GroupSpec:
    mode: str
    indices: Optional[Union[int, List[Union[int, str]]]] = None
    range: Optional[Union[str, List[int], Tuple[int, int]]] = None
    solute_end: Optional[int] = None
 
 
@dataclass(frozen=True)
class WallSpec:
    kind: str  # "U" or "L"
    kspring: float
    value: float
    power: int = 2
 
 
@dataclass(frozen=True)
class PlumedCVConfig:
    system: str
    buffer: float
    prefix: str
    method_dir: str
    bench_tag: str
    cv_dir: str
    replica_dir: str
    traj_name: str
    run_ids: List[int]
 
    cv_label: str
    # cv block
    dt_fs: float
    refdist: float
    nexp: int
    mexp: int
    norm: str
    # metadynamics block
    gausswidth: float
    metaheight_hartree: float
    metafreq_fs: float
    colvar_stride_fs: float
    grid_min: float
    grid_max: float
    grid_step: float
    group1: GroupSpec
    group2: GroupSpec
    wall: Optional[WallSpec] = None
 
 
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
 
 
def parse_indices(value: Union[int, List[Union[int, str]]]) -> List[int]:
    if isinstance(value, int):
        return [value]
    if not isinstance(value, list):
        raise RuntimeError(f"Bad indices spec: {value!r}")
    out: List[int] = []
    for item in value:
        if isinstance(item, int):
            out.append(item)
        elif isinstance(item, str) and "-" in item:
            a, b = parse_range(item)
            out.extend(range(a, b + 1))
        else:
            raise RuntimeError(f"Bad index item: {item!r}")
    return out
 
 
def select_group_indices(group: GroupSpec, symbols: List[str]) -> List[int]:
    """Return 1-based atom indices (PLUMED convention)."""
    mode = group.mode.lower()
    if mode == "indices":
        if group.indices is None:
            raise RuntimeError("group.indices required for mode=indices")
        return parse_indices(group.indices)
    if mode == "range":
        if group.range is None:
            raise RuntimeError("group.range required for mode=range")
        a, b = parse_range(group.range)
        return list(range(a, b + 1))
    if mode in ("all_water_h", "all_h"):
        if group.solute_end is None:
            raise RuntimeError("group.solute_end required for mode=all_water_H")
        start = int(group.solute_end) + 1
        out = [i for i, s in enumerate(symbols, start=1)
               if i >= start and s.upper() == "H"]
        if not out:
            raise RuntimeError("mode=all_water_H found no solvent hydrogens")
        if group.indices is None:
            return out
        extra = parse_indices(group.indices)
        seen, merged = set(), []
        for i in extra + out:
            if i not in seen:
                seen.add(i)
                merged.append(i)
        return merged
    raise RuntimeError(f"Unknown group.mode: {group.mode!r}")
 
 
def load_group(d: Dict[str, Any]) -> GroupSpec:
    return GroupSpec(
        mode=str(d["mode"]),
        indices=d.get("indices"),
        range=d.get("range"),
        solute_end=d.get("solute_end"),
    )
 
 
# --------------------------------------------------------------------------- #
# config / layout
# --------------------------------------------------------------------------- #
def load_config(yaml_path: Path) -> PlumedCVConfig:
    data = yaml.safe_load(yaml_path.read_text())
    cv = data["cv"]
    md = data["metadynamics"]
 
    wall_cfg = None
    if md.get("wall") is not None:
        w = md["wall"]
        kind = str(w["kind"]).strip().upper()
        if kind not in ("U", "L"):
            raise RuntimeError("wall.kind must be 'U' or 'L'")
        wall_cfg = WallSpec(kind=kind, kspring=float(w["kspring"]),
                            value=float(w["value"]), power=int(w.get("power", 2)))
 
    return PlumedCVConfig(
        system=data["system"],
        buffer=float(data["buffer"]),
        prefix=data.get("prefix", "solv"),
        method_dir=data.get("method_dir", "orb"),
        bench_tag=data["bench_tag"],
        cv_dir=data["cv_dir"],
        replica_dir=data.get("replica_dir", "equil"),
        traj_name=data.get("traj_name", "traj.xyz"),
        run_ids=parse_run_ids(data["run_ids"]),
        cv_label=data.get("cv_label", "c1"),
        dt_fs=float(cv["dt_fs"]),
        refdist=float(cv["refdist"]),
        nexp=int(cv["nexp"]),
        mexp=int(cv["mexp"]),
        norm=str(cv.get("norm", "AVERAGE")).upper(),
        gausswidth=float(md["gausswidth"]),
        metaheight_hartree=float(md["metaheight_hartree"]),
        metafreq_fs=float(md["metafreq_fs"]),
        colvar_stride_fs=float(md.get("colvar_stride_fs", md["metafreq_fs"])),
        grid_min=float(md["grid_min"]),
        grid_max=float(md["grid_max"]),
        grid_step=float(md["grid_step"]),
        group1=load_group(md["group1"]),
        group2=load_group(md["group2"]),
        wall=wall_cfg,
    )
 
 
class _orb_view:
    """Shim so we can reuse mlip.orb.orb_root_dir (expects .orb_dirname)."""
    def __init__(self, cfg: PlumedCVConfig):
        self.system = cfg.system
        self.buffer = cfg.buffer
        self.prefix = cfg.prefix
        self.orb_dirname = cfg.method_dir
 
 
def bench_dir(cfg: PlumedCVConfig, repo_root: Path) -> Path:
    return orb_root_dir(_orb_view(cfg), repo_root) / cfg.bench_tag
 
 
# --------------------------------------------------------------------------- #
# plumed.dat emission
# --------------------------------------------------------------------------- #
def build_plumed_lines(cfg: PlumedCVConfig, g1: List[int], g2: List[int]) -> List[str]:
    pace = max(1, round(cfg.metafreq_fs / cfg.dt_fs))
    print_stride = max(1, round(cfg.colvar_stride_fs / cfg.dt_fs))
    height_ev = cfg.metaheight_hartree * HARTREE_TO_EV
    nbins = max(1, round((cfg.grid_max - cfg.grid_min) / cfg.grid_step))
    npairs = len(g1) * len(g2)
 
    ga = ",".join(map(str, g1))
    gb = ",".join(map(str, g2))
 
    lines = ["UNITS LENGTH=A TIME=fs ENERGY=eV"]
    lines.append(
        f"{cfg.cv_label}: COORDINATION GROUPA={ga} GROUPB={gb} "
        f"R_0={cfg.refdist} NN={cfg.nexp} MM={cfg.mexp} D_0=0.0"
    )
 
    arg = cfg.cv_label
    if cfg.norm == "AVERAGE" and npairs > 1:
        arg = "cv"
        lines.append(f"cv: CUSTOM ARG={cfg.cv_label} FUNC=x/{npairs} PERIODIC=NO")
 
    if cfg.wall is not None:
        action = "UPPER_WALLS" if cfg.wall.kind == "U" else "LOWER_WALLS"
        lines.append(
            f"wall: {action} ARG={arg} AT={cfg.wall.value} "
            f"KAPPA={cfg.wall.kspring} EXP={cfg.wall.power}"
        )
 
    lines.append(
        f"metad: METAD ARG={arg} PACE={pace} HEIGHT={height_ev:.8f} "
        f"SIGMA={cfg.gausswidth} GRID_MIN={cfg.grid_min} GRID_MAX={cfg.grid_max} "
        f"GRID_BIN={nbins} FILE=HILLS"
    )
    lines.append(f"PRINT ARG={arg},metad.bias STRIDE={print_stride} FILE=COLVAR")
    return lines
 
 
def read_symbols(equil_dir: Path, traj_name: str) -> List[str]:
    traj = equil_dir / traj_name
    if not traj.exists() or traj.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty trajectory: {traj}")
    return list(read(traj, index=-1).get_chemical_symbols())
 
 
def _fmt_sel(tag: str, idx: List[int], symbols: List[str]) -> str:
    pairs = ", ".join(f"{i}:{symbols[i - 1]}" for i in idx)
    return f"{tag} {idx} -> [{pairs}]"
 
 
def run_cv(yaml_path: Path) -> None:
    cfg = load_config(yaml_path)
    repo_root = find_repo_root(yaml_path)
    bench = bench_dir(cfg, repo_root)
    yaml_text = yaml_path.read_text()
 
    for run_id in cfg.run_ids:
        run_root = bench / f"run-{run_id}"
        equil_dir = run_root / cfg.replica_dir
        cv_dir = run_root / cfg.cv_dir
        if (cv_dir / "plumed.dat").exists():
            print(f"SKIP: {cv_dir}/plumed.dat exists; not touching")
            continue
 
        symbols = read_symbols(equil_dir, cfg.traj_name)
        natoms = len(symbols)
        g1 = select_group_indices(cfg.group1, symbols)
        g2 = select_group_indices(cfg.group2, symbols)
        for i in g1 + g2:
            if i < 1 or i > natoms:
                raise RuntimeError(f"Atom index {i} out of range 1..{natoms} (run {run_id})")
 
        # confirmation: print 1-based index -> element so the user can eyeball it
        print(f"run {run_id}: {_fmt_sel('GROUPA', g1, symbols)} | "
              f"{_fmt_sel('GROUPB', g2, symbols)}")
 
        cv_dir.mkdir(parents=True, exist_ok=True)
        lines = build_plumed_lines(cfg, g1, g2)
        (cv_dir / "plumed.dat").write_text("\n".join(lines) + "\n")
        (cv_dir / "cv_spec.yaml").write_text(yaml_text)
        print(f"OK: wrote {cv_dir}/plumed.dat")
 
 
def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python3 plumed.py configs/SYS/plumed/nb.yaml", file=sys.stderr)
        raise SystemExit(2)
    run_cv(Path(sys.argv[1]).resolve())
 
 
if __name__ == "__main__":
    main()
