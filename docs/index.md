# mlip

Minimal, MLIP‑first workflows for molecular simulation. The goal is to keep setup and execution predictable while surfacing the ML interatomic potential (MLIP) runs front‑and‑center.

## Quickstart

```bash
conda env create -f environment.yml
conda activate mlip
pip install -e .
mlip --help
```

## Typical Flow

1. `prep` — build the solvated system with AmberTools/tleap
2. `mdequil` — equilibrate with classical force field (optional)
3. `salt` — adjust ions / hydroxide (optional)
4. `mlip-prep` — prepare MLIP runs
5. `mlip-submit` — submit MLIP runs to Slurm

## Design Notes

- MLIP runs are the center of the workflow; prep exists to make them reliable.
- Each command is YAML‑driven and writes outputs under `systems/<system>/<prefix>_<buffer>/`.
- `mlip` is independent from `qmmd` and shares no runtime linkage.
