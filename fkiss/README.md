# fkiss

`fkiss` (**F**ortran **KISS** — "Keep It Simple, Stupid") is a lightweight,
regex-based parser and dependency-analysis library for ABINIT's Fortran
source tree. It does not do compiler-level semantic analysis: it scans
`.F90`/`.finc` files with regular expressions to identify modules,
subroutines/functions, datatypes, interfaces, and `use`/CPP-driven
dependencies between them, and builds an in-memory model of the whole
`src/` tree from that.

It is a library, not a standalone tool — it has no `__main__` entry point of
its own. It is driven by **[`abisrc.py`](../abisrc.py)**, the CLI script at
the top of the ABINIT repository.

## Relation to `abisrc.py`

`abisrc.py` is the thin CLI front-end; `fkiss` is where all the actual
parsing and analysis logic lives:

- `fkiss.parser.FortranKissParser` / `FortranFile` — parse a single Fortran
  source file into modules, procedures, datatypes and their dependencies.
- `fkiss.project.AbinitProject` — parses (or unpickles a cached parse of)
  the whole `src/` tree into a dependency graph (DAG) and exposes the
  higher-level operations `abisrc.py`'s subcommands wrap: looking up a
  file/module/procedure/datatype, validating the tree, generating the
  `config/specs/*.conf` build-system files, drawing dependency graphs
  (via `graphviz`), and computing which files need to be recompiled after
  an API change (`touch`).
- `fkiss.regex` — the shared library of compiled regexes both `parser.py`
  and `project.py` build on.
- `fkiss.tools` — small shared utilities (`lazy_property`, dataframe/table
  printing, a `NotebookWriter` mixin used to open a project in Jupyter).
- `fkiss.termcolor` — ANSI terminal coloring, independent copy of the same
  helper used elsewhere in the ABINIT tooling (e.g.
  `tests/pymods/termcolor.py`).
- `fkiss.list_cpp_options` / `fkiss.mkrobodoc_dirs` / `fkiss.check_linalg_calls`
  — standalone analyses invoked from a few of `abisrc.py`'s subcommands
  (`cpp`, and ad hoc checks) or run directly.
- `fkiss.viewer` — an optional `panel`-based dashboard (`abisrc.py panel`)
  for browsing the project interactively in a browser.

`AbinitProject` parses the whole tree once and pickles the result (see
`AbinitProject.pickle_dump()`/`pickle_load()`) next to the source tree, so
that repeated `abisrc.py` invocations don't have to re-parse thousands of
Fortran files every time; `needs_reload()` compares file mtimes against the
pickle to decide whether a fresh parse is required.

## Typical usage

`fkiss` itself is only ever used as a library — always go through
`abisrc.py` from the top of the ABINIT source tree:

```bash
# Parse/refresh the dependency graph (auto-triggered by other commands too;
# rarely called directly except to force a rebuild with -r).
./abisrc.py parse .

# Inspect a file, module, datatype, interface, or public procedure.
./abisrc.py print src/41_geometry/m_crystal.F90
./abisrc.py print m_crystal
./abisrc.py print crystal_t
./abisrc.py print crystal_init

# Draw a dependency graph (requires graphviz + the python graphviz package).
./abisrc.py graph src/41_geometry/m_crystal.F90

# Developer workflows.
./abisrc.py makemake     # (Re)generate config/specs/*.conf build-system files.
./abisrc.py validate     # Validate the source tree (undeclared deps, etc.).
./abisrc.py touch        # Touch files (and their parents) after an API change,
                          # so `make` knows to recompile them.
./abisrc.py stats        # Pandas DataFrame with per-file/per-dir statistics.
./abisrc.py orphans      # Public procedures/modules that nobody depends on.
./abisrc.py cpp          # List CPP options used across the tree (fkiss.list_cpp_options).
./abisrc.py panel        # Interactive dashboard (fkiss.viewer, requires `panel`).
```

Run `./abisrc.py -h` (or `./abisrc.py <command> -h`) for the full, current
list of subcommands and options — `abisrc.py`'s own `--help` epilog is the
source of truth and can drift from this file.

## Tests

Unit tests live alongside the modules they cover (`fkiss/test_*.py`) and run
with `pytest` from the repository root, e.g.:

```bash
pytest fkiss
```
