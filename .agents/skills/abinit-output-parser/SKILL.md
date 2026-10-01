---
name: abinit-output-parser
description: Inspect an ABINIT text output or supported NetCDF result with AbiPy and summarize energies, convergence, forces, or stress. Use for ABINIT results, not arbitrary logs.
---

# Abinit Output Parser

This skill extracts key physical quantities and calculation metrics from an Abinit `.abo`/`.out`
text output file or a netCDF output file (e.g. `*_GSR.nc`).

## Usage

```bash
python .agents/skills/abinit-output-parser/scripts/parse_out.py --file "abinit_scf.abo"
```

### Arguments

- `--file`: Path to the Abinit output file (`.abo`/`.out` text log, or `.nc`, e.g. `*_GSR.nc`).

### Returns

A JSON object containing the parsed properties. **What is available depends on the file type**:

- **Plain text `.abo`/`.out`**: `energy_eV` (final SCF energy of the last cycle in the file),
  `num_scf_iterations`, and `converged` (`false` if abipy's event parser found an
  `Scf`/`Nscf`/`RelaxConvergenceWarning`, with the warning text under
  `convergence_warnings`). Text output does **not** expose structured forces/stress through
  abipy's API — parse the corresponding `*_GSR.nc` file for those (produced automatically
  alongside a ground-state run).
- **netCDF (`GSR.nc` and similar)**: `energy_eV`, `cart_forces_eV_per_ang`,
  `cart_stress_tensor_GPa` (the latter is `None`/omitted for NSCF runs, since stress is
  only meaningful for a self-consistent run).

A text output can contain several datasets.
The current helper reports only the last non-empty ground-state SCF cycle, so state this limitation
when presenting the result.
Missing JSON keys mean that AbiPy did not expose the quantity; they do not mean zero.
Do not infer scientific correctness solely from the absence of convergence warnings.

If a `.abo` file records an SCF cycle with zero iterations (e.g. a degenerate
`nstep=0`/single-point entry), that cycle is skipped rather than raising an error; if no
cycle with iterations is found, `energy_eV`/`num_scf_iterations` are simply omitted from
the result rather than crashing.

### Prerequisites

- `abipy` Python library must be installed.
