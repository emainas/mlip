#!/usr/bin/env python3
"""ParmEd-backed worker for :mod:`mlip.qmmm`.

This file is invoked with Amber's bundled Python so the main MLIP environment
does not need ParmEd as a Python dependency.
"""

import argparse
import json
from pathlib import Path

import parmed as pmd
from parmed.constants import PrmtopPointers
from parmed.topologyobjects import AtomType
from parmed.tools.actions import add12_6_4, addLJType
from parmed.tools.argumentlist import ArgumentList


def _read_ints(path: Path) -> list[int]:
    return [int(line) for line in path.read_text().splitlines() if line.strip()]


def _remove_periodic_box(parm: pmd.Structure) -> None:
    parm.box = None
    if hasattr(parm, "parm_data"):
        pointers = parm.parm_data.get("POINTERS")
        if pointers is not None:
            pointers[PrmtopPointers.IFBOX] = 0
        for flag in ("SOLVENT_POINTERS", "ATOMS_PER_MOLECULE", "BOX_DIMENSIONS"):
            parm.parm_data.pop(flag, None)
            if hasattr(parm, "flag_list") and flag in parm.flag_list:
                parm.flag_list.remove(flag)
            if hasattr(parm, "formats"):
                parm.formats.pop(flag, None)
            if hasattr(parm, "parm_comments"):
                parm.parm_comments.pop(flag, None)


def crop(args: argparse.Namespace) -> None:
    source = pmd.load_file(args.topology, xyz=args.coordinates)
    keep_residues = set(_read_ints(Path(args.keep_residues)))
    qm_water_residues = (
        set(_read_ints(Path(args.tip3p_qm_water_residues)))
        if args.tip3p_qm_water_residues
        else set()
    )

    # Selected OPC waters become three-site TIP3P waters before entering QM.
    # Giving the oxygen a distinct atom type preserves OPC parameters for the
    # surrounding MM solvent while using TIP3P LJ parameters at the boundary.
    type_number = len({atom.type for atom in source.atoms}) + 1
    tip3p_oxygen_type = AtomType("QO", type_number, 16.0, atomic_number=8)
    tip3p_oxygen_type.set_lj_params(0.1520, 1.7683)
    selected_water_atoms = set()
    removed_extra_points = 0
    for residue in source.residues:
        residue_number = residue.idx + 1
        if residue_number not in qm_water_residues:
            continue
        residue.name = "WAT"
        names = {atom.name for atom in residue.atoms}
        if not {"O", "H1", "H2"}.issubset(names):
            raise ValueError(
                f"QM water residue {residue_number} lacks O/H1/H2 atoms: {sorted(names)}"
            )
        for atom in residue.atoms:
            if atom.name == "O":
                atom.type = "QO"
                atom.atom_type = tip3p_oxygen_type
                atom.charge = -0.834
                selected_water_atoms.add(atom.idx)
            elif atom.name in {"H1", "H2"}:
                atom.charge = 0.417
                selected_water_atoms.add(atom.idx)
            else:
                removed_extra_points += 1

    keep_atoms = [
        atom.idx
        for residue in source.residues
        if residue.idx + 1 in keep_residues
        for atom in residue.atoms
        if not (residue.idx + 1 in qm_water_residues and atom.name not in {"O", "H1", "H2"})
    ]
    qm_water_indices = [
        new_index
        for new_index, source_index in enumerate(keep_atoms)
        if source_index in selected_water_atoms
    ]
    carved = source[keep_atoms]
    _remove_periodic_box(carved)
    if qm_water_residues:
        addLJType(
            carved,
            ArgumentList("@%QO radius 1.7683 epsilon 0.1520"),
        ).execute()
    if args.ion_1264_mask:
        action = add12_6_4(
            carved,
            ArgumentList(
                f"{args.ion_1264_mask} watermodel {args.ion_1264_watermodel}"
            ),
        )
        action.execute()
    carved.save(args.output_topology, overwrite=True)
    carved.save(args.output_coordinates, overwrite=True)
    print(
        json.dumps(
            {
                "atoms": len(carved.atoms),
                "residues": len(carved.residues),
                "total_partial_charge": sum(atom.charge for atom in carved.atoms),
                "qm_water_residues_source": sorted(qm_water_residues),
                "qm_water_atoms_zero_based": qm_water_indices,
                "removed_qm_water_extra_points": removed_extra_points,
            }
        )
    )


def inspect(args: argparse.Namespace) -> None:
    parm = pmd.load_file(args.topology, xyz=args.coordinates)
    qm = set(_read_ints(Path(args.qm_indices)))
    links = []
    for bond in parm.bonds:
        in1 = bond.atom1.idx in qm
        in2 = bond.atom2.idx in qm
        if in1 == in2:
            continue
        qatom, matom = (bond.atom1, bond.atom2) if in1 else (bond.atom2, bond.atom1)
        links.append(
            {
                "qm_index_zero_based": qatom.idx,
                "qm_atom": qatom.name,
                "qm_residue": qatom.residue.name,
                "qm_residue_number": qatom.residue.idx + 1,
                "mm_index_zero_based": matom.idx,
                "mm_atom": matom.name,
                "mm_residue": matom.residue.name,
                "mm_residue_number": matom.residue.idx + 1,
            }
        )
    print(
        json.dumps(
            {
                "atoms": len(parm.atoms),
                "residues": len(parm.residues),
                "qm_atoms": len(qm),
                "link_bonds": links,
            }
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    crop_parser = sub.add_parser("crop")
    crop_parser.add_argument("--topology", required=True)
    crop_parser.add_argument("--coordinates", required=True)
    crop_parser.add_argument("--keep-residues", required=True)
    crop_parser.add_argument("--output-topology", required=True)
    crop_parser.add_argument("--output-coordinates", required=True)
    crop_parser.add_argument("--ion-1264-mask")
    crop_parser.add_argument("--ion-1264-watermodel", default="OPC")
    crop_parser.add_argument("--tip3p-qm-water-residues")
    crop_parser.set_defaults(func=crop)

    inspect_parser = sub.add_parser("inspect")
    inspect_parser.add_argument("--topology", required=True)
    inspect_parser.add_argument("--coordinates", required=True)
    inspect_parser.add_argument("--qm-indices", required=True)
    inspect_parser.set_defaults(func=inspect)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
