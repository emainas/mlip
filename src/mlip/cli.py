import argparse
from pathlib import Path

from mlip.prep import run_prep
from mlip.mdequil import run_mdequil
from mlip.salt import run_salt
from mlip.rdf import run_rdf
from mlip.dihedral import run_dihedral
from mlip.rmsd import run_rmsd
from mlip.radgyr import run_radgyr
from mlip.hbond import run_hbond
from mlip.orb import run_orb_prep, run_orb_submit
from mlip.plumed import run_cv
from mlip.orb_meta import run_orb_meta_prep, run_orb_meta_submit


def main():
    p = argparse.ArgumentParser(prog="mlip")
    sub = p.add_subparsers(dest="cmd", required=True)

    prep = sub.add_parser("prep", help="Prepare solvated system with tleap")
    prep.add_argument("yaml", type=Path)

    mdequil = sub.add_parser("mdequil", help="Write MD equil inputs and run (slurm if provided, else local)")
    mdequil.add_argument("yaml", type=Path)

    salt = sub.add_parser("salt", help="Delete the counterion and turn the furthest water into a hydroxide")
    salt.add_argument("yaml", type=Path)


    orb_prep = sub.add_parser("orb-prep", help="Write ORB equil inputs/scripts (no submit)")
    orb_prep.add_argument("yaml", type=Path)

    orb_submit = sub.add_parser("orb-submit", help="Submit ORB equil jobs for runs matching the config")
    orb_submit.add_argument("yaml", type=Path)

    rdf = sub.add_parser("rdf", help="Compute RDF with cpptraj for selected runs")
    rdf.add_argument("yaml", type=Path)

    dihedral = sub.add_parser("dihedral", help="Compute dihedral time series with cpptraj for selected runs")
    dihedral.add_argument("yaml", type=Path)

    rmsd = sub.add_parser("rmsd", help="Compute RMSD time series with cpptraj for selected runs")
    rmsd.add_argument("yaml", type=Path)

    radgyr = sub.add_parser("radgyr", help="Compute radius of gyration time series with cpptraj for selected runs")
    radgyr.add_argument("yaml", type=Path)

    hbond = sub.add_parser("hbond", help="Compute hydrogen bond time series and lifetimes with cpptraj")
    hbond.add_argument("yaml", type=Path)

    plumed = sub.add_parser("plumed", help="Write plumed.dat for selected runs")
    plumed.add_argument("yaml", type=Path)

    orb_meta_prep = sub.add_parser("orb-meta-prep", help="Write ORB+PLUMED meta inputs/scripts (no submit)")
    orb_meta_prep.add_argument("yaml", type=Path)

    orb_meta_submit = sub.add_parser("orb-meta-submit", help="Submit ORB+PLUMED meta runs for matching configs")
    orb_meta_submit.add_argument("yaml", type=Path)

    args = p.parse_args()

    if args.cmd == "prep":
        run_prep(args.yaml)
    elif args.cmd == "mdequil":
        run_mdequil(args.yaml)
    elif args.cmd == "salt":
        run_salt(args.yaml)
    elif args.cmd == "orb-prep":
        run_orb_prep(args.yaml)
    elif args.cmd == "orb-submit":
        run_orb_submit(args.yaml)
    elif args.cmd == "rdf":
        run_rdf(args.yaml)
    elif args.cmd == "dihedral":
        run_dihedral(args.yaml)
    elif args.cmd == "rmsd":
        run_rmsd(args.yaml)
    elif args.cmd == "radgyr":
        run_radgyr(args.yaml)
    elif args.cmd == "hbond":
        run_hbond(args.yaml)
    elif args.cmd == "plumed":
        run_cv(args.yaml)
    elif args.cmd == "orb-meta-prep":
        run_orb_meta_prep(args.yaml)
    elif args.cmd == "orb-meta-submit":
        run_orb_meta_submit(args.yaml)


if __name__ == "__main__":
    main()
