# QM/MM snapshot preparation

Prepare non-periodic, spherical Amber snapshots for TeraChem's built-in
OpenMM QM/MM driver.

```bash
PYTHONPATH=src python -m mlip.cli qmmm configs/osteocalcin/qmmm/qmmm.yaml
```

With the package installed, the equivalent command is `mlip qmmm
configs/osteocalcin/qmmm/qmmm.yaml`.

The osteocalcin YAML treats the 10 ns production NetCDF as the mother run. It
requires all 1000 saved frames to be present, then selects frames 1, 112, 223,
334, 445, 556, 667, 778, 889, and 1000. The check prevents a partially written
trajectory from changing the meaning of `stop: last`.

Each snapshot is imaged around the five calcium ions and converted to a
non-periodic 28 Å sphere containing the full dimer and all neutralizing sodium
ions. The charge -4 singlet QM region contains all five calcium ions, the side
chains of GLA 6, 10, 13, 44, 48, and 51, ASP 19 and 57, and every water oxygen
within 2.8 Å of a calcium. Coordinating OPC waters are converted to three-site
TIP3P waters before their O/H atoms enter QM; their EP sites are removed. This
gives 116-119 QM atoms in the sampled frames. Eight CA-CB boundaries are
validated for TeraChem's automatic link atoms.

Each frame directory receives a carved `system.parm7`, `frame.rst7`, zero-based
`region.qm`, `tc.in`, and a validation manifest. The `terachem` YAML mapping is
written directly into `tc.in`, so the calculation type is selected through the
input rather than a separate workflow.

For a single-point energy with pure spherical basis functions:

```yaml
terachem:
  run: energy
  sphericalbasis: "yes"
```

For QM/MM molecular dynamics, retain `run: md` and place the MD controls under
`extra`:

```yaml
terachem:
  run: md
  seed_base: 20261026
  extra:
    nstep: 1000
    timestep: 0.5
    velocities: random
    tinit: 300.0
    thermostat: rescale
    mdbc: spherical
    md_density: 1.0
```

The packaged osteocalcin configuration retains the 500 fs MD setup for later
use. The currently submitted calculation uses the single-point form above.
TeraChem 25.03's embedded Amber reader
only reads the Lennard-Jones A and B tables, however, so it ignores the C4 term
during QM/MM. All five calcium ions are in the QM region, which reduces but does
not eliminate the consequence for Ca interactions with MM water. After
reviewing the manifests, run the generated calculations from the output
directory with `./run_all.sh`.
