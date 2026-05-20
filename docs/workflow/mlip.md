# mlip

Prepare and submit MLIP (ORB) runs.

## Commands

```bash
mlip mlip-prep configs/<molecule>/mlip/mlip.yaml
mlip mlip-submit configs/<molecule>/mlip/mlip.yaml
```

## Outputs

- `run.sh`, `slurm.sh` (if configured)
- `main.py` runner
- `spec.yaml`
