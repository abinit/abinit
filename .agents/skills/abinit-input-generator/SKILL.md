---
name: abinit-input-generator
description: Generate a starting ABINIT .abi input from a real structure and local pseudopotentials with AbiPy. Use for input scaffolding, not converged production settings.
---

# Abinit Input Generator

This skill generates a valid Abinit `.abi` input file for a given structure and calculation type.

The result is a starting point for review, not a scientifically validated production input.

## Usage

```bash
python .agents/skills/abinit-input-generator/scripts/generate.py \
    --structure Si.cif --calc-type scf --pseudos-dir /path/to/pseudos --output Si_scf.abi
```

### Arguments

- `--structure`: Path to a structure file (`.cif`, POSCAR, etc.). Bare chemical formulas
  (e.g. `"Si"`) are **not** supported: there is no reliable, offline way to reconstruct a
  real crystal structure (lattice parameters, space group) from a formula alone, and
  earlier versions of this skill silently fabricated a bogus single-atom cubic cell in
  that case, which produced physically meaningless input files without any warning.
- `--calc-type`: `scf`, `relax`, `nscf`, or `bands`. `nscf` and `bands` both produce a
  2-dataset SCF+NSCF `MultiDataset` via abipy's `ebands_input` factory — abipy has no
  bare, single-dataset NSCF factory that doesn't start from a prior SCF density, so there
  is no meaningful difference between the two here.
- `--pseudos-dir`: Directory containing the pseudopotential file for every element in the
  structure (e.g. a PseudoDojo table you already downloaded locally). Required — there is
  no automatic online pseudopotential lookup.
- `--output`: Path to write the generated `.abi` file.
  Do not overwrite an existing file unless the user explicitly authorizes replacement.

### Prerequisites

- `abipy` (and `pymatgen`) must be installed.
- The pseudopotential files for the structure's elements must already exist locally.

After generation, inspect the structure, pseudopotential mapping, cutoff, k-point sampling,
occupations, convergence criteria, and requested outputs.
Validate the final input with the intended ABINIT executable.
