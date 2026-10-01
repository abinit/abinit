---
name: abinit-add-fortran-file
description: Add or register a Fortran source file in ABINIT's generated Autotools and CMake metadata. Use for new files under src/ or shared/*/src, not ordinary edits to registered files.
---

# Add Fortran File to Build System

This skill automates the registration of a new Fortran source file in the Abinit build system.
When you create a new `.F90` file inside a directory (e.g. `src/78_eph/m_my_module.F90`),
this skill updates `abinit.src` and `CMakeLists.txt` in that directory to include the new file.

Work from the repository root and read `config/README.md` before changing generated build metadata.
Use the helper only when inspection confirms that both local files use the supported `sources = [...]`
and `add_library(...)` layouts.

## Usage

```bash
python .agents/skills/abinit-add-fortran-file/scripts/add_fortran_file.py --file src/78_eph/m_my_module.F90
```

### Arguments

- `--file`: The relative path to an existing `.F90` file under `src/`, `shared/common/src/`,
  or `shared/libpaw/src/`.

### What it does

1. Validates that the file and directory exist.
2. Reads `abinit.src` in that directory.
3. Finds the `sources = [...]` array and inserts the file name alphabetically.
4. Reads `CMakeLists.txt` in that directory.
5. Finds the `add_library(...)` source list and inserts the file name alphabetically.

After running the helper, inspect `git diff` immediately.
Preserve alphabetical ordering and all existing comments.
If the directory has multiple CMake targets or an unfamiliar layout, edit deliberately instead of
forcing the helper.

Run `./config/scripts/makemake`, reconfigure the build directory when necessary, compile the affected
target, and run relevant tests.
Treat any helper error as a failed operation and inspect for partial changes before retrying.
