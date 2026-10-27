# dihedral

Dihedral time‑series analysis with `cpptraj`.

```bash
mlip dihedral configs/<molecule>/analysis/dihedral.yaml
```

## Plotting

```bash
python plotting/plot_dihedral_timeseries.py \
  --runs-path systems/<system>/<prefix>_<buffer>/orb/<bench_tag> \
  --labels double5,single5,single10,double10,single15,double15 \
  --window 50 \
  --out reports/dihedrals.png
```
