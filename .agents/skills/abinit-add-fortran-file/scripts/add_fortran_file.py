#!/usr/bin/env python3
"""
Registers a new Fortran file in the abinit build system.
Updates abinit.src and CMakeLists.txt in the file's directory.
"""

import argparse
import os
import sys


def insert_alphabetically(lines, start_idx, end_idx, new_line, extract_key):
    """
    Inserts a new_line alphabetically into lines between start_idx and end_idx.
    extract_key is a function to extract the string to compare from a line.
    """
    new_key = extract_key(new_line)
    insert_idx = start_idx
    for i in range(start_idx, end_idx):
        key = extract_key(lines[i])
        if key > new_key:
            break
        insert_idx = i + 1

    lines.insert(insert_idx, new_line)
    return lines

def update_abinit_src(directory, filename):
    abinit_src = os.path.join(directory, "abinit.src")
    if not os.path.exists(abinit_src):
        print(f"Error: {abinit_src} not found.")
        return False

    with open(abinit_src) as f:
        lines = f.readlines()

    start_idx = -1
    end_idx = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("sources = ["):
            start_idx = i + 1
        elif start_idx != -1 and line.strip() == "]":
            end_idx = i
            break

    if start_idx == -1 or end_idx == -1:
        print(f"Error: Could not find 'sources = [' block in {abinit_src}.")
        return False

    extract_key = lambda x: x.strip().strip('",')

    # Check if already exists
    for i in range(start_idx, end_idx):
        if extract_key(lines[i]) == filename:
            print(f"Info: {filename} already exists in {abinit_src}.")
            return True

    new_line = f'"{filename}",\n'

    insert_alphabetically(lines, start_idx, end_idx, new_line, extract_key)

    with open(abinit_src, "w") as f:
        f.writelines(lines)
    print(f"✅ Added {filename} to {abinit_src}")
    return True

def update_cmake(directory, filename):
    cmake_file = os.path.join(directory, "CMakeLists.txt")
    if not os.path.exists(cmake_file):
        print(f"Error: {cmake_file} not found.")
        return False

    with open(cmake_file) as f:
        lines = f.readlines()

    start_idx = -1
    end_idx = -1
    for i, line in enumerate(lines):
        if "add_library(" in line:
            start_idx = i + 1
        elif start_idx != -1 and line.strip() == ")":
            end_idx = i
            break

    if start_idx == -1 or end_idx == -1:
        print(f"Error: Could not find 'add_library(' block in {cmake_file}.")
        return False

    extract_key = lambda x: x.strip()

    # Check if already exists
    for i in range(start_idx, end_idx):
        if extract_key(lines[i]) == filename:
            print(f"Info: {filename} already exists in {cmake_file}.")
            return True

    new_line = f"  {filename}\n"

    insert_alphabetically(lines, start_idx, end_idx, new_line, extract_key)

    with open(cmake_file, "w") as f:
        f.writelines(lines)
    print(f"✅ Added {filename} to {cmake_file}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Register a new Fortran file in the build system.")
    parser.add_argument("--file", type=str, required=True, help="Path to the new Fortran file (e.g. src/78_eph/m_my_mod.F90)")
    args = parser.parse_args()

    repository_root = os.path.realpath(os.getcwd())
    file_path = os.path.realpath(args.file)
    allowed_roots = [
        os.path.join(repository_root, "src"),
        os.path.join(repository_root, "shared", "common", "src"),
        os.path.join(repository_root, "shared", "libpaw", "src"),
    ]
    if not os.path.isfile(file_path) or not file_path.endswith(".F90"):
        print(f"Error: Expected an existing .F90 file, got {args.file}")
        return 1
    if not any(os.path.commonpath([file_path, root]) == root for root in allowed_roots):
        print(f"Error: {args.file} is outside an ABINIT source directory.")
        return 1

    directory = os.path.dirname(file_path)
    filename = os.path.basename(file_path)

    if not directory:
        print("Error: Please provide the relative path including the directory, e.g. src/78_eph/m_mod.F90")
        return 1

    print(f"Registering {filename} in {directory}...")

    metadata_paths = [os.path.join(directory, name) for name in ("abinit.src", "CMakeLists.txt")]
    if not all(os.path.isfile(path) for path in metadata_paths):
        print(f"Error: Expected abinit.src and CMakeLists.txt in {directory}.")
        return 1

    originals = {}
    for path in metadata_paths:
        with open(path, encoding="utf-8") as fh:
            originals[path] = fh.read()

    if not (update_abinit_src(directory, filename) and update_cmake(directory, filename)):
        for path, content in originals.items():
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(content)
        print("Error: Registration failed; restored both metadata files.")
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
