#!/usr/bin/env python3
"""
Parse an Abinit output file using abipy.
"""

import argparse
import json
import os
import sys

try:
    from abipy.abilab import abiopen
    from abipy.abio.outputs import AbinitOutputFile
    from abipy.flowtk import events as abievents
    from pymatgen.core.units import Energy
except ImportError:
    print(json.dumps({"error": "abipy is not installed."}))
    sys.exit(1)

# Event classes that mean "this run did not converge", across the calculation
# types abiopen's AbinitOutputFile.events can report on.
CONVERGENCE_WARNING_CLASSES = (
    abievents.ScfConvergenceWarning,
    abievents.NscfConvergenceWarning,
    abievents.RelaxConvergenceWarning,
)


def parse_abo(abifile):
    """
    Extract what is actually available from a plain-text .abo/.out file.
    Note: unlike netCDF output (e.g. GSR.nc), plain text output does not expose
    structured forces/stress via abipy's API, so those are not returned here.
    """
    data = {}

    cycles = abifile.get_all_gs_scf_cycles()
    # Some cycles (e.g. nstep=0 single-point runs) record zero iterations, so
    # "Etot(hartree)" is empty and last_etotal would raise -- skip those.
    cycles = [c for c in cycles if c.num_iterations > 0]
    if cycles:
        last_cycle = cycles[-1]
        data["energy_eV"] = float(Energy(last_cycle.last_etotal, "Ha").to("eV"))
        data["num_scf_iterations"] = last_cycle.num_iterations

    warnings = [e for e in abifile.events if isinstance(e, CONVERGENCE_WARNING_CLASSES)]
    data["converged"] = len(warnings) == 0
    if warnings:
        data["convergence_warnings"] = [str(w) for w in warnings]

    return data


def parse_ncfile(abifile):
    """Extract energy/forces/stress from a netCDF output file (e.g. GSR.nc)."""
    data = {}
    if hasattr(abifile, "energy"):
        data["energy_eV"] = float(abifile.energy)
    if hasattr(abifile, "cart_forces") and abifile.cart_forces is not None:
        data["cart_forces_eV_per_ang"] = abifile.cart_forces.tolist()
    if hasattr(abifile, "cart_stress_tensor") and abifile.cart_stress_tensor is not None:
        data["cart_stress_tensor_GPa"] = abifile.cart_stress_tensor.tolist()
    return data


def main():
    parser = argparse.ArgumentParser(description="Parse Abinit output file.")
    parser.add_argument("--file", type=str, required=True, help="Path to Abinit output (.abo/.out or .nc) file.")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(json.dumps({"error": f"File not found: {args.file}"}))
        sys.exit(1)

    try:
        with abiopen(args.file) as abifile:
            if isinstance(abifile, AbinitOutputFile):
                data = parse_abo(abifile)
            else:
                data = parse_ncfile(abifile)
            print(json.dumps({"status": "success", "data": data}))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)

if __name__ == "__main__":
    main()
