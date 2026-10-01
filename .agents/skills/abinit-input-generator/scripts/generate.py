#!/usr/bin/env python3
"""
Generate an Abinit input file using abipy.
"""

import argparse
import glob
import json
import os
import sys

try:
    import abipy.abio.factories as abif
    from abipy.core.structure import Structure
    from pymatgen.io.abinit.pseudos import PseudoTable
except ImportError:
    print(json.dumps({"error": "abipy is not installed. Please install it to use this skill."}))
    sys.exit(1)


def load_pseudos(pseudos_dir):
    """Build a PseudoTable from every pseudopotential file in pseudos_dir."""
    paths = sorted(
        p for p in glob.glob(os.path.join(pseudos_dir, "*"))
        if os.path.isfile(p) and not p.endswith((".json", ".md", ".txt"))
    )
    if not paths:
        raise ValueError(f"No pseudopotential files found in {pseudos_dir}")
    return PseudoTable.as_table(paths)


def main():
    parser = argparse.ArgumentParser(description="Generate Abinit input file.")
    parser.add_argument("--structure", type=str, required=True,
                         help="Path to a structure file (e.g. .cif, POSCAR). "
                              "Bare chemical formulas are not supported: there is no reliable, "
                              "offline way to look up a real crystal structure from a formula alone.")
    parser.add_argument("--calc-type", type=str, choices=["scf", "nscf", "relax", "bands"],
                         default="scf",
                         help="Calculation type. 'nscf' and 'bands' both produce a 2-dataset "
                              "SCF+NSCF MultiDataset (abipy has no bare, single-dataset NSCF "
                              "factory that doesn't start from a prior SCF density).")
    parser.add_argument("--pseudos-dir", type=str, required=True,
                         help="Directory containing the pseudopotential files for every element "
                              "in the structure (e.g. a PseudoDojo table you already downloaded).")
    parser.add_argument("--output", type=str, required=True, help="Path to output .abi file.")

    args = parser.parse_args()

    if os.path.exists(args.output):
        print(json.dumps({"error": f"Output file already exists: {args.output}"}))
        sys.exit(1)

    if not os.path.exists(args.structure):
        print(json.dumps({
            "error": f"Structure file not found: {args.structure}. "
                     "Pass a real structure file (.cif, POSCAR, etc.); "
                     "bare formulas like 'Si' cannot be turned into a correct structure offline."
        }))
        sys.exit(1)
    structure = Structure.from_file(args.structure)

    try:
        pseudos = load_pseudos(args.pseudos_dir)
    except Exception as e:
        print(json.dumps({"error": f"Failed to load pseudopotentials: {e!s}"}))
        sys.exit(1)

    try:
        if args.calc_type == "scf":
            inp = abif.scf_input(structure, pseudos)
        elif args.calc_type == "relax":
            inp = abif.ion_ioncell_relax_input(structure, pseudos)
        elif args.calc_type in ("nscf", "bands"):
            inp = abif.ebands_input(structure, pseudos)
        else:
            raise ValueError(f"Unhandled calc-type: {args.calc_type}")

        with open(args.output, "w") as f:
            f.write(str(inp))

        print(json.dumps({"status": "success", "output_file": args.output}))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)

if __name__ == "__main__":
    main()
