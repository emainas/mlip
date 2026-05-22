#!/usr/bin/env python3

import os
import sys
import yaml
import shutil
import subprocess 
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PrepConfig:
    """Class for spefifying the configuration of a solvation prep with tleap (Ambertools)."""
    start_from: str
    name: str
    leaprc_mol: str
    leaprc_sol: str
    frcmod: Path
    mol2: Path
    water_model: str
    buffer: float
    input_dir: Path
    xyz: Optional[Path] = None
    template_mol2: Optional[Path] = None
    template_frcmod: Optional[Path] = None
    charge: Optional[int] = None
    charge_method: Optional[str] = None
    resname: Optional[str] = None
    counterion_num: int = 0
    frcmod_ion: Optional[str] = None
    counterion: Optional[str] = None
    prefix: Optional[str] = None

def load_config(yaml_path: Path) -> PrepConfig:
    
    data = yaml.safe_load(yaml_path.read_text())
    
    cfg = PrepConfig(
        start_from=data.get("start_from", "mol2"),
        name=data["name"],
        leaprc_mol=data["leaprc_mol"],
        leaprc_sol=data["leaprc_sol"],
        frcmod_ion=data["frcmod_ion"],
        frcmod=Path(data["frcmod"]),
        mol2=Path(data["mol2"]),
        water_model=data.get("water_model"),
        buffer=float(data.get("buffer")),
        counterion=data.get("counterion"),
        counterion_num=int(data.get("counterion_num")),
        input_dir=Path(data["input_dir"]),
        prefix=data.get("prefix", "solv"),
        xyz=Path(data["xyz"]).expanduser() if data.get("xyz") else None,
        template_mol2=Path(data["template_mol2"]).expanduser() if data.get("template_mol2") else None,
        template_frcmod=Path(data["template_frcmod"]).expanduser() if data.get("template_frcmod") else None,
        charge=int(data["charge"]) if data.get("charge") is not None else None,
        charge_method=data.get("charge_method"),
        resname=data.get("resname"),
    )

    return cfg

def read_xyz(xyz_path: Path) -> tuple[list[str], list[tuple[float, float, float]]]:
    with xyz_path.open("r") as f:
        first = f.readline()
        try:
            natoms = int(first.strip())
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ header in {xyz_path}") from exc
        f.readline()  # comment
        symbols = []
        coords = []
        for _ in range(natoms):
            line = f.readline()
            if not line:
                break
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ line in {xyz_path}: {line.rstrip()}")
            symbols.append(parts[0])
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    if len(coords) != natoms:
        raise ValueError(f"Expected {natoms} atoms in {xyz_path}, got {len(coords)}")
    return symbols, coords

def generate_atom_names(symbols: list[str]) -> list[str]:
    atom_types: dict[str, int] = {}
    names = []
    for sym in symbols:
        count = atom_types.get(sym, 0) + 1
        atom_types[sym] = count
        token = sym + "%0" + str(3 - len(sym)) + "d"
        names.append(token % count)
    return names

def read_mol2_template(
    mol2_path: Path,
) -> tuple[list[dict], list[tuple[int, int, str]], list[str], str, str]:
    atoms = []
    bonds = []
    substructure_lines = []
    resname = None
    mol_name = None
    in_atoms = False
    in_bonds = False
    in_sub = False
    with mol2_path.open("r") as f:
        for line in f:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol_name = f.readline().strip()
                continue
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                in_bonds = False
                in_sub = False
                continue
            if line.startswith("@<TRIPOS>BOND"):
                in_atoms = False
                in_bonds = True
                in_sub = False
                continue
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                in_atoms = False
                in_bonds = False
                in_sub = True
                continue
            if line.startswith("@<TRIPOS>") and not line.startswith("@<TRIPOS>ATOM") and not line.startswith("@<TRIPOS>BOND") and not line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                in_atoms = False
                in_bonds = False
                in_sub = False
            if in_atoms:
                parts = line.split()
                if len(parts) < 9:
                    continue
                atom_id = int(parts[0])
                atoms.append(
                    {
                        "id": atom_id,
                        "name": parts[1],
                        "type": parts[5],
                        "resid": int(parts[6]),
                        "resname": parts[7],
                        "charge": parts[8],
                    }
                )
                if resname is None:
                    resname = parts[7]
            if in_bonds:
                parts = line.split()
                if len(parts) < 4:
                    continue
                bonds.append((int(parts[1]), int(parts[2]), parts[3]))
            if in_sub:
                if line.strip():
                    substructure_lines.append(line.rstrip("\n"))
    if not atoms:
        raise ValueError(f"No atoms found in template mol2: {mol2_path}")
    return atoms, bonds, substructure_lines, resname or "MOL", mol_name or "MOL"

def write_mol2_from_template(
    dest: Path,
    mol_name: str,
    atoms: list[dict],
    bonds: list[tuple[int, int, str]],
    substructure_lines: list[str],
    symbols: list[str],
    coords: list[tuple[float, float, float]],
) -> None:
    if len(atoms) != len(coords):
        raise ValueError("Template mol2 atom count does not match XYZ atom count")
    with dest.open("w") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{mol_name}\n")
        f.write(f"{len(atoms):5d}{len(bonds):6d}     1     0     0\n")
        f.write("SMALL\n")
        f.write("USER_CHARGES\n\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for atom, sym, (x, y, z) in zip(atoms, symbols, coords):
            f.write(
                "%7d %-8s %10.4f %10.4f %10.4f %-6s %4d %-8s %10s\n"
                % (
                    atom["id"],
                    atom["name"],
                    x,
                    y,
                    z,
                    atom["type"],
                    atom["resid"],
                    atom["resname"],
                    atom["charge"],
                )
            )
        f.write("@<TRIPOS>BOND\n")
        for i, (a1, a2, btype) in enumerate(bonds, start=1):
            f.write(f"{i:6d} {a1:5d} {a2:5d} {btype}\n")
        if substructure_lines:
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            for line in substructure_lines:
                f.write(f"{line}\n")

def write_pdb_simple(
    dest: Path,
    atom_names: list[str],
    resname: str,
    symbols: list[str],
    coords: list[tuple[float, float, float]],
) -> None:
    if len(atom_names) != len(coords):
        raise ValueError("Atom name count does not match XYZ atom count")
    with dest.open("w") as f:
        f.write(f"COMPND   {resname}\n")
        for i, (name, sym, (x, y, z)) in enumerate(zip(atom_names, symbols, coords), start=1):
            f.write(
                "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                % (i, name, resname, " ", 1, x, y, z, 1.00, 0.00, sym)
            )

def read_mol2_atoms(mol2_path: Path) -> list[dict]:
    atoms = []
    in_atoms = False
    with mol2_path.open("r") as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>") and in_atoms:
                break
            if in_atoms:
                parts = line.split()
                if len(parts) < 9:
                    continue
                atoms.append(
                    {
                        "id": int(parts[0]),
                        "name": parts[1],
                        "x": float(parts[2]),
                        "y": float(parts[3]),
                        "z": float(parts[4]),
                        "type": parts[5],
                        "resid": int(parts[6]),
                        "resname": parts[7],
                        "charge": parts[8],
                    }
                )
    if not atoms:
        raise ValueError(f"No atoms found in mol2: {mol2_path}")
    return atoms

def require_ambertools_tool(name: str) -> Path:
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        raise RuntimeError("AMBERHOME is not set.")

    exe = shutil.which(name)
    if not exe:
        raise RuntimeError(f"{name} not found in PATH. Did you load AmberTools?")

    exe_path = Path(exe).resolve()
    expected = Path(amberhome).resolve() / "bin" / name
    if exe_path != expected:
        raise RuntimeError(
            f"{name} mismatch:\n"
            f"  PATH {name}: {exe_path}\n"
            f"  AMBERHOME: {amberhome}\n"
            f"  expected: {expected}\n"
            f"Fix by re-loading the correct Amber module."
        )

    return exe_path

def run_antechamber(cfg: PrepConfig, prep_dir: Path) -> None:
    if not cfg.xyz:
        raise RuntimeError("xyz is required for start_from: xyz")
    if cfg.charge is None or not cfg.charge_method:
        raise RuntimeError("charge and charge_method are required for start_from: xyz")

    xyz_path = cfg.xyz.resolve()
    symbols, coords = read_xyz(xyz_path)

    antechamber = require_ambertools_tool("antechamber")
    antechamber_log = prep_dir / "antechamber.out"

    if cfg.template_mol2:
        template_mol2 = cfg.template_mol2.resolve()
        atoms, bonds, sub_lines, template_resname, mol_name = read_mol2_template(template_mol2)
        resname = cfg.resname or template_resname
        for atom in atoms:
            atom["resname"] = resname

        merged_mol2 = prep_dir / (xyz_path.stem + "_template.mol2")
        write_mol2_from_template(merged_mol2, mol_name, atoms, bonds, sub_lines, symbols, coords)
        antechamber_out = prep_dir / f"{cfg.mol2.stem}_antechamber.mol2"
        with antechamber_log.open("w") as f:
            subprocess.run(
                [
                    str(antechamber),
                    "-i",
                    merged_mol2.name,
                    "-fi",
                    "mol2",
                    "-o",
                    antechamber_out.name,
                    "-fo",
                    "mol2",
                    "-c",
                    cfg.charge_method,
                    "-nc",
                    str(cfg.charge),
                    "-rn",
                    resname,
                ],
                cwd=prep_dir,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )

        ac_atoms = read_mol2_atoms(antechamber_out)
        if len(ac_atoms) != len(atoms):
            raise RuntimeError(
                f"Antechamber atom count {len(ac_atoms)} does not match template {len(atoms)}"
            )

        final_mol2 = prep_dir / cfg.mol2.name
        with final_mol2.open("w") as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{mol_name}\n")
            f.write(f"{len(atoms):5d}{len(bonds):6d}     1     0     0\n")
            f.write("SMALL\n")
            f.write("USER_CHARGES\n\n\n")
            f.write("@<TRIPOS>ATOM\n")
            for tmpl, ac, (x, y, z) in zip(atoms, ac_atoms, coords):
                f.write(
                    "%7d %-8s %10.4f %10.4f %10.4f %-6s %4d %-8s %10s\n"
                    % (
                        tmpl["id"],
                        tmpl["name"],
                        x,
                        y,
                        z,
                        ac["type"],
                        tmpl["resid"],
                        resname,
                        ac["charge"],
                    )
                )
            f.write("@<TRIPOS>BOND\n")
            for i, (a1, a2, btype) in enumerate(bonds, start=1):
                f.write(f"{i:6d} {a1:5d} {a2:5d} {btype}\n")
            if sub_lines:
                f.write("@<TRIPOS>SUBSTRUCTURE\n")
                for line in sub_lines:
                    f.write(f"{line}\n")
    else:
        resname = cfg.resname or "MOL"
        atom_names = generate_atom_names(symbols)
        pdb_path = prep_dir / (xyz_path.stem + ".pdb")
        write_pdb_simple(pdb_path, atom_names, resname, symbols, coords)
        with antechamber_log.open("w") as f:
            subprocess.run(
                [
                    str(antechamber),
                    "-i",
                    pdb_path.name,
                    "-fi",
                    "pdb",
                    "-o",
                    cfg.mol2.name,
                    "-fo",
                    "mol2",
                    "-c",
                    cfg.charge_method,
                    "-nc",
                    str(cfg.charge),
                    "-rn",
                    resname,
                ],
                cwd=prep_dir,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )

    if cfg.template_frcmod and cfg.template_frcmod.exists():
        shutil.copy2(cfg.template_frcmod, prep_dir / cfg.frcmod.name)
        return

    parmchk2 = require_ambertools_tool("parmchk2")
    parmchk_log = prep_dir / "parmchk2.out"
    with parmchk_log.open("w") as f:
        subprocess.run(
            [
                str(parmchk2),
                "-i",
                cfg.mol2.name,
                "-o",
                cfg.frcmod.name,
                "-f",
                "mol2",
                "-s",
                "2",
                "-a",
                "Y",
            ],
            cwd=prep_dir,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
        )

def write_tleap_in(cfg: PrepConfig) -> Path:

    base_dir = cfg.input_dir.parent
    prep_dir = base_dir / f"{cfg.prefix}_{cfg.buffer:.1f}" / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)

    frcmod_name = cfg.frcmod.name
    mol2_name = cfg.mol2.name

    frcmod_src = cfg.input_dir / frcmod_name
    mol2_src = cfg.input_dir / mol2_name
    
    # Existence of input checks
    if cfg.start_from == "mol2":
        if not frcmod_src.exists():
            raise FileNotFoundError(f"Missing in input_dir: {frcmod_src}")
        if not mol2_src.exists():
            raise FileNotFoundError(f"Missing in input_dir: {mol2_src}")
    
    # Copy inputs to prep dir to run tleap in
    if cfg.start_from == "mol2":
        shutil.copy2(frcmod_src, prep_dir / frcmod_name)
        shutil.copy2(mol2_src, prep_dir / mol2_name)
    elif cfg.start_from == "xyz":
        if cfg.template_mol2 and cfg.template_frcmod and cfg.template_frcmod.exists():
            symbols, coords = read_xyz(cfg.xyz.resolve())
            atoms, bonds, sub_lines, template_resname, mol_name = read_mol2_template(cfg.template_mol2.resolve())
            resname = cfg.resname or template_resname
            for atom in atoms:
                atom["resname"] = resname
            write_mol2_from_template(prep_dir / mol2_name, mol_name, atoms, bonds, sub_lines, symbols, coords)
            shutil.copy2(cfg.template_frcmod, prep_dir / frcmod_name)
        else:
            run_antechamber(cfg, prep_dir)
    else:
        raise ValueError(f"Unsupported start_from: {cfg.start_from}")

    tleap_in = prep_dir / "tleap.in"

    addions_block = ""
    if cfg.counterion_num > 0:
        addions_block = f'addions sys {cfg.counterion} {cfg.counterion_num}\n'

    text = f"""\
# Auto-generated by mlip prep.py
# System: {cfg.name}

source {cfg.leaprc_mol}
source {cfg.leaprc_sol}

# Ions (Joung-Cheatham for TIP3P)
loadamberparams {cfg.frcmod_ion}

# Molecule params
loadamberparams {frcmod_name}
sys = loadmol2 {mol2_name}

# Solvate (buffer in Angstrom)
solvatebox sys {cfg.water_model} {cfg.buffer}

# Counterions
{addions_block.rstrip()}

saveamberparm sys {cfg.prefix}.parm7 {cfg.prefix}.rst7
quit
"""
    tleap_in.write_text(text)
    return tleap_in


def require_amber() -> Path:
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        raise RuntimeError("AMBERHOME is not set.")

    tleap = shutil.which("tleap")
    if not tleap:
        raise RuntimeError("tleap not found in PATH. Did you load Amber?")

    tleap_path = Path(tleap).resolve()

    expected = Path(amberhome).resolve() / "bin" / "tleap"
    if tleap_path != expected:
        raise RuntimeError(
            f"tleap mismatch:\n"
            f"  PATH tleap: {tleap_path}\n"
            f"  AMBERHOME: {amberhome}\n"
            f"  expected: {expected}\n"
            f"Fix by re-loading the correct Amber module."
        )

    return tleap_path

def run_tleap(cfg: PrepConfig, tleap_in: Path) -> None: 
    
    prep_dir = tleap_in.parent
    out_path = prep_dir / "tleap.out"

    with out_path.open("w") as f:
        tleap = require_amber()
        subprocess.run(
            [str(tleap), "-f", tleap_in.name],
            cwd=prep_dir,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
        )

def run_prep(yaml_path):
    cfg = load_config(yaml_path)
    tleap_in = write_tleap_in(cfg)
    run_tleap(cfg, tleap_in)

    prep_dir = tleap_in.parent

    # Freeze config for reproducibility (in the prep dir)
    (prep_dir / "spec.yaml").write_text(yaml_path.read_text())

    parm7 = prep_dir / f"{cfg.prefix}.parm7"
    rst7 = prep_dir / f"{cfg.prefix}.rst7"
    if not parm7.exists() or parm7.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty output: {parm7}")
    if not rst7.exists() or rst7.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty output: {rst7}")

    print(f"OK: wrote {parm7.name} and {rst7.name} in {prep_dir}")

def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python3 prep.py prep.yaml", file=sys.stderr)
        raise SystemExit(2)

    yaml_path = Path(sys.argv[1]).resolve()
    cfg = load_config(yaml_path)

    tleap_in = write_tleap_in(cfg)
    run_tleap(cfg, tleap_in)

    prep_dir = tleap_in.parent

    # Freeze config for reproducibility (in the prep dir)
    (prep_dir / "spec.yaml").write_text(yaml_path.read_text())

    parm7 = prep_dir / f"{cfg.prefix}.parm7"
    rst7 = prep_dir / f"{cfg.prefix}.rst7"
    if not parm7.exists() or parm7.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty output: {parm7}")
    if not rst7.exists() or rst7.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty output: {rst7}")

    print(f"OK: wrote {parm7.name} and {rst7.name} in {prep_dir}")


if __name__ == "__main__":
    main()
