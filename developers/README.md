# ABINIT developer and maintainer utilities

This directory contains scripts used for source maintenance, release preparation, test-reference management, debugging, and one-off migrations.
These are repository-maintenance tools, not runtime components of ABINIT.

The directory includes a mixture of active tools and historical utilities.
Do not assume that every executable file represents a supported workflow.
Read a script before running it, check its assumptions against the current source tree, and review all resulting changes with Git.

General development instructions are available in the [developer overview](../doc/developers/overview_development.md) and the [developer how-to](../doc/developers/developers_howto.md).
The test infrastructure is documented in [tests/README.md](../tests/README.md), and the generated build-system machinery is described in [config/README.md](../config/README.md).

## Directory map

| Path | Purpose |
|---|---|
| `maintainers/` | Repository-wide checks and transformations, release chores, and reference-file maintenance. |
| `various/` | Developer convenience scripts, debugging aids, and historical conversion tools. |
| `pltdiff.py` | Interactive plotting of numeric test output against reference data. |

None of these directories is a Python package or a stable public API.
Scripts generally expect to be run from a particular directory and often derive paths from the current working directory.

## Safety first

Several tools modify many tracked files and have no dry-run mode.
Before using one of them:

1. Start from a clean or deliberately isolated Git worktree.
2. Read the complete script and identify its hard-coded patterns, paths, and temporary filenames.
3. Run it first on the smallest practical set of files.
4. Inspect `git status` and `git diff` immediately.
5. Run `./config/scripts/makemake` if source files or build metadata changed.
6. Run the relevant tests and pre-commit hooks before keeping the result.

Do not run repository-rewriting scripts concurrently.
Some older shell utilities use shared temporary names such as `tmp.file`, and an interrupted run can leave partial changes behind.

The strongest caution applies to:

- `change.sh`, `suppress.sh`, and `change_year.sh`, which replace files in place;
- `replace_write6.py` and `suppress_robodoc_section.py`, which rewrite Fortran sources;
- `check_datatypes.py --autofix`, which edits `use` statements;
- `remove_unused_var.sh`, which parses a compiler log and edits declarations;
- `update_refs.sh`, which copies generated test results over reference files and removes comparison reports.

Git is the recovery mechanism for these operations.
Never run them on uncommitted work unless those changes are understood and backed up.

## Maintainer tools

### `abirules.pl`

`maintainers/abirules.pl` is the principal tool in this directory that is integrated with the generated Makefiles.
It checks and can rewrite ABINIT Fortran source layout, RoboDOC headers, declarations, `use` statements, preprocessor sections, and selected coding conventions.

Its command-line form is:

```sh
perl developers/maintainers/abirules.pl [options]
```

Important options include:

| Option | Effect |
|---|---|
| `-d PATH` | Process one source file or directory instead of the default source-tree selection. |
| `-p` | Preserve originals and write files with an `.abirules` suffix. |
| `-r` | Do not reorder declarations. |
| `-l` | Reorder declarations line by line. |
| `-v` | Print verbose diagnostics. |

Prefer `-p` or a narrowly scoped `-d` run while evaluating a change.
The default behavior can replace source files.

The build-system generators in `config/scripts/make-makefiles-inter` create `abirules` targets in intermediate Makefiles.
The selected directories come from the build-system configuration and are processed by this Perl script.
Therefore, changes to the Make interface belong in the generator under `config/scripts`, not only in a generated `Makefile.am`.

Run the generated target from a configured build tree with:

```sh
make abirules
```

Review `tmp-abirules.log` and all modified files afterward.
Out-of-source behavior should be checked carefully because parts of the generated recipe still use paths relative to source directories.

### `change_year.sh`

`maintainers/change_year.sh` performs the annual ABINIT copyright update.
It is referenced by the release checklist in `doc/maintainers/check_list.txt`.

The replacement years, filename patterns, exceptional paths, and final permission changes are hard-coded in the script.
At the start of a new year, read and update the script before applying it from the repository root.
Its header contains the intended multi-command workflow and the searches needed to find missed copyright notices.

This is not a generic `old-year new-year` command.
It removes and recreates files while processing them and applies broad `chmod` commands at the end, so the resulting content and modes both require review.

### `update_refs.sh`

`maintainers/update_refs.sh` bulk-updates test reference files from `TestBot_MPI*` result directories.
It is intended for an authorized reference machine after a deliberate decision to accept new numerical results.

Run it from the top-level `tests` directory, as stated in its header.
The suite names and output globs are hard-coded and must be reviewed whenever suites or reference artifacts change.
After copying outputs, the script removes matching `.fldiff` reports so that remaining reports reveal outputs not covered by its mapping.

Updating reference files is not a routine test-repair operation.
Inspect scientific differences, compiler and library changes, test status, and the final Git diff before accepting new references.

### Source checks and transformations

| Script | Purpose | Status and cautions |
|---|---|---|
| `check-missing-headers.py` | Reports coverage of RoboDOC module and routine headers in a source tree. | Read-only; defaults to `src`. Its parser is heuristic. |
| `check_datatypes.py` | Examines selected Fortran datatype modules and can suggest narrower `use` statements. | Uses hard-coded module/type inventories dating from an older source layout. Validate the inventory before use; `--autofix` rewrites files. |
| `replace_write6.py` | Replaces literal standard-output/error unit numbers with ABINIT unit variables. | Repository-wide source rewrite designed for a specific migration. Use only after auditing current matches. |
| `suppress_robodoc_section.py` | Removes a hard-coded RoboDOC section, currently `PARENTS`, from source files. | Destructive migration utility with no command-line section selector or dry run. |
| `change.sh` | Applies a hard-coded `sed` substitution to files given on the command line. | Despite its historical description, it does not accept search and replacement expressions and is currently effectively a template/no-op substitution. Edit and review it before use. |
| `suppress.sh` | Removes lines matching a hard-coded pattern, currently `src2tex`. | Destructive migration utility using a shared temporary file. |
| `fix-cpp-options.py` | Renames a fixed set of historical preprocessor symbols. | Legacy Python 2 code (`file()` and list-like dictionary keys); it is not usable unchanged with the supported Python environment. |
| `suppress_windows_ctrl-M.csh` | Removes Windows carriage returns from selected files. | Legacy C-shell cleanup tool; prefer a targeted modern line-ending check when possible. |

These scripts are not invoked by the normal compilation or test suite unless an explicit Make target or maintainer calls them.

## Developer convenience tools

### `pltdiff.py`

`pltdiff.py` plots columns from numeric output and reference files with NumPy and Matplotlib.
With one output path, it tries to locate the corresponding `tests/SUITE/Refs` file from the build-tree layout.
With two paths, it compares those files directly.

Examples:

```sh
python developers/pltdiff.py path/to/test-output
python developers/pltdiff.py reference-file output-file
```

Input must be whitespace-separated numeric tabular data with an initial x-axis column.
Blank lines and lines beginning with `#` are ignored.
The tool opens an interactive Matplotlib window and is not suitable for headless regression testing without additional backend configuration.

### Source scaffolding

`various/mkmodule.sh` and `various/mkroutine.sh` create Fortran skeletons in the current directory:

```sh
developers/various/mkmodule.sh m_new_module new_routine
developers/various/mkroutine.sh new_routine
```

The generated text contains deliberate `FIXME` markers and incomplete declarations that must be replaced.
The templates also reflect historical RoboDOC and source-layout conventions, so compare the result with a recently added neighboring source file.

After creating a source file, register it in both build frontends:

- update the authoritative Autotools source/dependency configuration and run `./config/scripts/makemake`;
- update the appropriate `CMakeLists.txt` source list and target dependencies.

Follow the source-registration instructions in the [developer how-to](../doc/developers/developers_howto.md).

### Parallel run and debugging helpers

| Script | Purpose | Limitations |
|---|---|---|
| `various/paradev.py` | Runs an executable with several MPI process counts and opens HTML diffs between outputs. | Assumes the historical `files` standard-input interface and uses `os.system`; its comparison is textual rather than `fldiff`-based. |
| `various/paragdb.py` | Starts one `xterm` and `gdb` process per MPI rank. | Requires `mpirun`, X11, `xterm`, and `gdb`; it is a crude helper, not a parallel debugger. |

Treat both as local interactive aids.
They are not portable automation interfaces and are not part of the test farm.

### Historical conversion and cleanup tools

| Script | Purpose | Current caveat |
|---|---|---|
| `various/fixed_to_free` | Converts fixed-form Fortran to free form. | Requires the external `FT77to90` program, which is not included in this directory despite the script's historical comment. |
| `various/remove_unused_var.sh` | Parses GNU Fortran unused-variable warnings from `make.log` and edits declarations. | Highly dependent on old compiler output and source formatting; it also leaves `.bak` files. Audit before any use. |

These utilities are retained as historical or specialized tools rather than recommended development workflows.

## Relationship with checks and tests

The normal contributor checks live primarily under `abichecks`, the Python test directories, and `.pre-commit-config.yaml`.
Use those maintained interfaces before reaching for an ad hoc script in this directory:

```sh
pre-commit run --all-files
pytest
```

For Fortran changes, run the relevant ABINIT test suites in addition to formatting and static checks.
`abirules.pl` can complement these checks but does not replace compilation or numerical testing.

The scripts here do not share a common test harness.
When modifying one, add a non-destructive mode and focused tests where practical, and make path assumptions explicit through command-line arguments rather than the current working directory.

## Build-system and distribution interfaces

The `developers` directory is not part of ABINIT's compilation or installation.
At present it is also not registered in `config/specs/buildsys.conf`, so ordinary `make dist` does not obtain it through a generated `config/dist/auto-developers.lst`.

There are nevertheless several legacy or maintainer-facing references:

- generated intermediate Makefiles call `developers/maintainers/abirules.pl` from their `abirules` target;
- `config/makefiles/top.am` has legacy binary-package rules that copy the `developers` directory;
- `dist-lite` removes `developers` from its temporary tree if it is present;
- the release checklist refers to `change_year.sh` and `update_refs.sh`.

Two generated `parents` targets still refer to `developers/maintainers/parents.pl`, but that script is no longer present.
Those targets are stale and should not be advertised or used until the generator is corrected.
Likewise, legacy binary-package recipes refer to maintainer Makefile templates that are no longer present in this directory.

If these utilities are intended to be shipped in the standard source archive, register `developers` as a detached data directory in `config/specs/buildsys.conf`, run `makemake`, and verify the resulting archive with `make dist`.
If they are intentionally development-tree-only, the existing exclusion should be documented in release policy and stale packaging references should be removed.

Do not edit generated `Makefile.am` or `Makefile.in` files to change these relationships.
Update `config/specs`, `config/makefiles`, or the generator under `config/scripts`, then run:

```sh
./config/scripts/makemake
```

See [config/README.md](../config/README.md) for the complete generated-file and distribution workflow.

## Improving or retiring a utility

Before extending an old script, decide whether its operation belongs in an existing maintained interface such as `abichecks`, `abisrc.py`, the test runner, or a pre-commit hook.
A maintained replacement should normally provide:

- `--help` with explicit working-directory and dependency requirements;
- a read-only check or `--dry-run` mode;
- explicit input paths and narrowly scoped output paths;
- atomic or recoverable writes;
- meaningful exit status and diagnostics;
- focused tests for parsing and transformations;
- Python 3 compatibility when implemented in Python.

When removing a script, search the build generators, Make fragments, release checklists, documentation, and CI configuration for references.
Run `makemake` if a generated Make interface changes, then execute the relevant checks and inspect `make dist` or legacy packaging targets as appropriate.
