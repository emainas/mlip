# rdf

Radial distribution function analysis with `cpptraj`.

```bash
mlip rdf configs/<molecule>/analysis/rdf.yaml
```

## Plotting

```bash
python plotting/rdf_hist.py \
  --root systems/<system>/<prefix>_<buffer>/orb/<bench_tag> \
  --label rdf_1_HD_WAT_O \
  --out reports/rdf_hist.png
```

```bash
python plotting/rdf_bv_hydrogens.py \
  --root systems/<system>/<prefix>_<buffer>/orb/<bench_tag> \
  --out reports/rdf_bv_hydrogens.png
```
