# mlip

Lightweight workflows for ML interatomic potential (MLIP) simulations, with Amber-based prep and analysis utilities.

---

## Installation

```bash
git clone <your-mlip-repo-url>
cd mlip
conda env create -f environment.yml
conda activate mlip
pip install -e .
mlip --help
```

---

## prep — System Preparation

Generates a solvated system using `tleap`.

```bash
mlip prep configs/<molecule>/prep/prep.yaml
```

---

## mdequil — MD Equilibration

NVT and NPT equilibration with classical force field using Amber's sander module.

```bash
mlip mdequil configs/<molecule>/mdequil/mdequil.yaml
```

---

## salt — Post-equilibration System Adjustment

Deletes tleap's counterion and mutates the furthest water molecule into hydroxide.

```bash
mlip salt configs/<molecule>/salt/salt.yaml
```

---

## density — Solute/Box Volume

Computes solute and total box volume from salt outputs.

```bash
mlip density configs/<molecule>/density/density.yaml
```

---

## mlip-prep — Prepare MLIP (ORB) Runs

Writes ORB run inputs/scripts (no submission).

```bash
mlip mlip-prep configs/<molecule>/mlip/mlip.yaml
```

---

## mlip-submit — Submit MLIP (ORB) Runs

Submits Slurm jobs for run directories that match the provided config.

```bash
mlip mlip-submit configs/<molecule>/mlip/mlip.yaml
```

---

## Analysis (cpptraj)

- `rdf`: radial distribution functions
- `dihedral`: dihedral time series
- `rmsd`: RMSD time series
- `radgyr`: radius of gyration time series
- `hbond`: hydrogen bond time series and lifetimes

```bash
mlip rdf configs/<molecule>/analysis/rdf.yaml
mlip dihedral configs/<molecule>/analysis/dihedral.yaml
mlip rmsd configs/<molecule>/analysis/rmsd.yaml
mlip radgyr configs/<molecule>/analysis/radgyr.yaml
mlip hbond configs/<molecule>/analysis/hbond.yaml
```
