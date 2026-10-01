#!/usr/bin/env python3
"""
Automates adding a new input variable to the Abinit Fortran codebase.
Based on doc/developers/developers_howto.md
"""

import argparse
import os
import re
import sys


def insert_after(file_path, pattern, text_to_insert):
    """Insert text after the first occurrence of a regex pattern in a file."""
    if not os.path.exists(file_path):
        return False

    with open(file_path) as f:
        content = f.read()

    match = re.search(pattern, content)
    if not match:
        return False

    insertion_point = match.end()
    # Add newline if text_to_insert doesn't start with one
    if not text_to_insert.startswith("\n"):
        text_to_insert = "\n" + text_to_insert

    new_content = content[:insertion_point] + text_to_insert + content[insertion_point:]

    with open(file_path, "w") as f:
        f.write(new_content)
    return True

def main():
    parser = argparse.ArgumentParser(description="Add a new Abinit input variable.")
    parser.add_argument("--name", type=str, required=True, help="Name of the variable (e.g. myvar)")
    parser.add_argument("--type", type=str, required=True, help="Fortran type (e.g. integer, real(dp))")
    parser.add_argument("--default", type=str, required=True, help="Default value")
    parser.add_argument("--description", type=str, required=True, help="Description for docs")

    args = parser.parse_args()

    var_name = args.name.lower()
    if not re.fullmatch(r"[a-z][a-z0-9_]*", var_name):
        print("Error: Variable name must be a valid lowercase-style identifier.")
        return 1
    if var_name[-1].isdigit():
        print("Error: Variable name cannot end with a digit.")
        return 1

    print(f"Adding variable '{var_name}' of type '{args.type}'...")

    # Check if we are in the root of the abinit repository
    if not os.path.exists("src/44_abitypes_defs/m_dtset.F90"):
        print("Error: Must be run from the root of the abinit repository.")
        return 1

    existing = []
    for path in (
        "src/44_abitypes_defs/m_dtset.F90",
        "src/57_iovars/m_invars1.F90",
        "abimkdocs/variables_abinit.py",
    ):
        with open(path, encoding="utf-8") as fh:
            if re.search(rf"\b{re.escape(var_name)}\b", fh.read(), flags=re.IGNORECASE):
                existing.append(path)
    if existing:
        print(f"Error: Variable {var_name!r} already appears in: {', '.join(existing)}")
        return 1

    success_files = []
    failed_files = []

    # 1. m_dtset.F90: dataset_type declaration + chkvars.
    # (There is no separate "dtset_copy" subroutine in the current codebase --
    # dtsets(:) inherit/duplicate defaults via m_invars1.F90's own dataset-0
    # propagation logic, not a per-field copy routine, so no step is needed here.)
    file1 = "src/44_abitypes_defs/m_dtset.F90"
    if insert_after(file1, r"type\s*,\s*public\s*::\s*dataset_type[\s\S]*?(?=end type dataset_type)", f"  {args.type} :: {var_name}\n"):
        success_files.append(f"{file1} (dataset_type)")
    else:
        failed_files.append(f"{file1} (dataset_type)")

    if insert_after(file1, r"subroutine\s+chkvars[\s\S]*?(?=end subroutine chkvars)", f"  ! TODO: Check {var_name}\n"):
         success_files.append(f"{file1} (chkvars)")
    else:
         failed_files.append(f"{file1} (chkvars)")

    # 2. m_invars1.F90: default value (indefo1) and read logic (invars1).
    file2 = "src/57_iovars/m_invars1.F90"
    if insert_after(file2, r"subroutine\s+indefo1[\s\S]*?(?=end subroutine indefo1)", f"  dtset%{var_name} = {args.default}\n"):
        success_files.append(f"{file2} (indefo1)")
    else:
        failed_files.append(f"{file2} (indefo1)")

    if insert_after(file2, r"subroutine\s+invars1\s*\([\s\S]*?(?=end subroutine invars1)", f"  ! TODO: call intagm(...) to read {var_name} from input\n"):
        success_files.append(f"{file2} (invars1)")
    else:
        failed_files.append(f"{file2} (invars1)")

    # 3. m_outvar_X.F90 (printing), split alphabetically like the real files.
    first_letter = var_name[0]
    if "a" <= first_letter <= "h":
        outvar_file, outvar_sub = "src/57_iovars/m_outvar_a_h.F90", "outvar_a_h"
    elif "i" <= first_letter <= "n":
        outvar_file, outvar_sub = "src/57_iovars/m_outvar_i_n.F90", "outvar_i_n"
    else:
        outvar_file, outvar_sub = "src/57_iovars/m_outvar_o_z.F90", "outvar_o_z"

    if insert_after(outvar_file, rf"subroutine\s+{outvar_sub}\s*\([\s\S]*?(?=end subroutine {outvar_sub})", f"  ! TODO: Print {var_name}\n"):
        success_files.append(f"{outvar_file} ({outvar_sub})")
    else:
        failed_files.append(f"{outvar_file} ({outvar_sub})")

    # 4. m_chkinp.F90: consistency checks.
    file4 = "src/57_iovars/m_chkinp.F90"
    if insert_after(file4, r"subroutine\s+chkinp\s*\([\s\S]*?(?=end subroutine chkinp)", f"  ! TODO: Check consistency of {var_name}\n"):
        success_files.append(file4)
    else:
        failed_files.append(file4)

    print("\n--- Summary ---")
    for f in success_files:
        print(f"✅ Modified: {f}")
    for f in failed_files:
        print(f"❌ Failed to modify: {f}")

    print("\nPlease manually review the changes using `git diff`.")
    if failed_files:
        print("One or more insertions failed; the working tree may contain partial scaffolding.")
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
