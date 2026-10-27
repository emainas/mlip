# Osteocalcin QM/MM sample

QM/MM combines two kinds of data: quantum-mechanical data for a chemically important inner region—here, atomic coordinates, elements, total charge, multiplicity, a list of 119 QM atoms containing five Ca²⁺ ions, coordinating carboxylates and first-shell waters, plus the electronic-structure method and basis—and molecular-mechanical data for the surrounding environment, namely atom types, bonded terms, partial charges, force-field parameters and coordinates for approximately 12,200 MM atoms; the calculation uses both regions together to obtain quantities such as the electronic energy and SCF convergence history. **Metadata** means “data about the data”: the contextual information needed to identify, interpret and reproduce a result, such as the source trajectory and frame, units, software version, QM/MM partition, force fields, charge and spin, hardware, run type, timestamps and final status. In this directory, `metadata.yaml` records that context in machine-readable form, while the remaining files provide a complete frame-1 input/output example.

## Included sample

- `tc.in` — TeraChem QM/MM single-point input.
- `system.parm7` — Amber topology, force-field parameters and MM charges.
- `frame.rst7` — coordinates for the cropped, nonperiodic snapshot.
- `region.qm` — zero-based indices of the 119 explicitly selected QM atoms.
- `system.pdb` — human-readable structure for inspection and visualization.
- `tc.out` — completed TeraChem output, including all 18 SCF iterations.
- `results.dat` — compact final-energy result written by TeraChem.
- `qm_region.png` — rendered view of the sampled QM region.
- `metadata.yaml` — provenance, model choices, dimensions, hardware and result summary.

With the TeraChem 25.03 environment and license available, the copied calculation can be rerun from this directory with:

```bash
terachem tc.in > rerun.out
```

