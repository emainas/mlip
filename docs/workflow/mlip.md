# mlip

Prepare and submit MLIP (ORB) runs.

## Commands

```bash
mlip orb-prep configs/<molecule>/orb/orb.yaml
mlip orb-submit configs/<molecule>/orb/orb.yaml
```

### Optional: PLUMED CV setup

```bash
mlip plumed-cv configs/<molecule>/orb-meta/plumed/plumed.yaml
```

## Outputs

- `run.sh`, `slurm.sh` (if configured)
- `main.py` runner
- `spec.yaml`

For `plumed-cv`:
- `plumed.dat`
- `cv_spec.yaml`
