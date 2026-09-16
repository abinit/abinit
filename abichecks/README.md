# ABINIT maintenance checks

The `abichecks` directory contains fast maintenance and source-consistency checks.
These checks complement the numerical and functional tests under `tests/`.

Most checks inspect source files, build-system specifications, generated files, or packaging rules rather than executing a scientific calculation.
They are organized into named series and can be launched through Make from a configured ABINIT build tree.

## Relationship with `make distcheck`

`make distcheck` is an Automake release-validation target.
It verifies that the source archive is complete and that the package can be built independently of the original Git checkout.

At a high level, `make distcheck` performs the following operations:

```text
create the source archive with make dist
                    |
                    v
             unpack the archive
                    |
                    v
       configure a separate build tree
                    |
                    v
          build and run make check
                    |
                    v
   install into a temporary staging area
                    |
                    v
       run installcheck and uninstall
                    |
                    v
 verify distclean and archive integrity
```

This process detects problems that an ordinary build from a Git checkout can hide.
For example, it detects required files missing from `EXTRA_DIST`, dependencies on files outside the archive, installation errors, and files left behind by `make distclean`.

The top-level `AM_DISTCHECK_CONFIGURE_FLAGS` variable supplies the configuration options used for the temporary build.
ABINIT disables several optional dependencies there so that release validation does not require every supported external library.

Run `distcheck` from a configured build directory:

```sh
source /Users/giantomassi/miniconda3/etc/profile.d/conda.sh
conda activate env3.14
cd _build
make distcheck
```

`make distcheck` is significantly slower than `make dist` because it configures, compiles, tests, installs, uninstalls, and cleans a fresh copy of the package.
It also needs enough free disk space for the archive, unpacked source, build products, and temporary installation.

### What `distcheck` does not imply

The custom targets documented below are not aliases for `make distcheck`.
Automake invokes the normal recursive `check` targets during `distcheck`, but it does not automatically invoke every project-specific target whose name starts with `tests_`.

In particular, run these explicitly when their coverage is required:

```sh
cd _build/abichecks
make tests_abirules
make tests_buildsys
```

`make distcheck` validates that the `abichecks` files are correctly distributed and that the standard check workflow succeeds.
The explicit commands above execute the three maintenance-check series.

The generated Makefile also provides `tests_libpaw`, and `run_checks.py` (the Python runner) dispatches it the same way as the other two series.

## Directory organization

The checks are divided into three series:

| Directory | Make target | Alias | Purpose |
|---|---|---|---|
| `abirules/` | `tests_abirules` | `Tabirules` | Enforce ABINIT source conventions and inspect compiler warnings. |
| `buildsys/` | `tests_buildsys` | `Tbuildsys` | Check build-system specifications and generated configuration. |
| `libpaw/` | `tests_libpaw` | `Tlibpaw` | Verify that LibPAW can be packaged and built independently. |

Each series contains a `tests.cnf` file.
The file assigns a numeric identifier to each command that belongs to the series.
Some series also contain `Input/` data and `Refs/` reference output used by the legacy comparison machinery.

The `scripts/` directory contains the implementations and the common runners.

## How Make targets are generated

The `abichecks/Makefile.am` file is generated and must not be edited directly.
The generation path is:

```text
config/specs/abichecks.conf
            |
            v
config/scripts/make-makefiles-abichecks
            |
            +---- reads config/dist/auto-abichecks.lst
            |
            v
     abichecks/Makefile.am
            |
            v
       Automake and configure
            |
            v
    _build/abichecks/Makefile
```

`config/specs/abichecks.conf` declares each series, its description, and optional Make aliases.
`config/scripts/make-makefiles-abichecks` converts these declarations into Make targets and adds the files listed in `config/dist/auto-abichecks.lst` to `EXTRA_DIST`.

The distribution list is produced by `config/scripts/make-file-lists` from the repository layout and the exclusions in `config/specs/junk.conf`.
The full process is normally driven by:

```sh
./config/scripts/makemake
```

When adding or removing an `abichecks` file, run `makemake` so that the generated distribution list, Makefile inputs, and Autotools outputs remain synchronized.

## How a check is executed

For example, this command:

```sh
cd _build/abichecks
make tests_buildsys
```

follows this path:

```text
tests_buildsys target in the configured Makefile
                        |
                        v
        . ./abichecks.env && $(PYTHON) scripts/run_checks.py
                        |
                        +---- reads abichecks/buildsys/tests.cnf
                        |
                        v
      dispatch each registered command or script
                        |
                        v
       write results under _build/abichecks
```

`abichecks.env` is generated from `abichecks.env.in` by `configure`.
It records the source directory, build directory, shell, and Python interpreter required by the runner.
This is why the Make targets should be run from the configured `abichecks` build directory rather than directly from the source directory.

By default, result directories use the hostname.
The `dirname` Make variable can select a different result-directory name.

## Selecting checks

Run a complete series with:

```sh
cd _build/abichecks
make tests_buildsys
```

Use `start` and `stop` to run an inclusive range of numeric check identifiers:

```sh
make tests_buildsys start=03 stop=12
```

To run one check, specify only `start`:

```sh
make tests_buildsys start=04
```

Use a custom output-directory name with:

```sh
make tests_buildsys start=04 dirname=my-check
```

The aliases are equivalent to the full target names:

```sh
make Tbuildsys
make Tabirules
make Tlibpaw
```

Run `make help` in `_build/abichecks` to display the concise target summary.

## The `tests.cnf` files

The `tests.cnf` files are the authoritative registration point for the legacy Make-driven series.
A typical entry is:

```text
04 statchk abichecks/scripts/check-unprotected-omp-target.py
```

The fields mean:

1. `04` is the identifier used by `start` and `stop`.
2. `statchk` tells the runner to execute a status-returning check.
3. The final field identifies the script to execute.

The `abirules/tests.cnf` series also uses commands such as `warnchk` and finishes with `make report` to compare the collected warning summaries with reference data.
Read the header of an existing `tests.cnf` before introducing a new command type or multi-step check.

## Script overview

### Common orchestration and comparison scripts

| Script | Role |
|---|---|
| `run_checks.py` | Pure-Python dispatcher, invoked directly by the Makefile recipe. Reads `tests.cnf`, prepares the `tmp-<host>_<os>_<date>` working directory, runs each registered check, and writes the `report` file plus per-case logs. Reuses `tests/pymods/fldiff.py`'s `Differ`/`Result` (an existing, tested floating-point-aware comparison engine) for golden-file checks instead of a separate diff tool. |
| `run-basic-tests.pl` | Legacy controller for built-in executable tests. |
| `run-basic-dotest.pl` | Executes one built-in test on behalf of `run-basic-tests.pl`. |
| `run-basic-tests.sh` | Shell implementation of the older built-in-test workflow. |

The basic-test scripts are retained for compatibility with older Make workflows and are unrelated to `run_checks.py`/`tests.cnf` above.
New numerical regression tests should normally use `tests/runtests.py` and the `TEST_INFO` format documented in `tests/README.md`.

### Source-rule checks

| Script | Role |
|---|---|
| `check_ascii.py` | Reports non-ASCII characters in supported source-file types. |
| `check_forbidden.py` | Detects forbidden source constructs such as explicit output units, direct allocation, and disallowed MPI communicators. |
| `check_inlined_macros.py` | Finds ABINIT memory macros embedded in unsupported inline expressions. |
| `warningschk.py` | Scans compiler-warning logs by warning category and applies the configured exclusions. |
| `abirules_tools.py` | Provides shared helpers for locating the repository and ABINIT source directories. |

These checks form the `abirules` series declared in `abirules/tests.cnf`.

### Build-system checks

| Script | Role |
|---|---|
| `check-config-h.py` | Verifies the standard guarded inclusion of `config.h` in preprocessed Fortran sources. |
| `check-unprotected-omp-target.py` | Detects OpenMP target regions that are not protected as required; it uses `gfortran -E` only for candidate files. |
| `check-build-config.py` | Checks consistency among build configuration keywords, environment definitions, and options. |
| `check-cpp-options.py` | Cross-checks C-preprocessor symbols used in source, headers, and Autoconf definitions. |
| `check-binaries-conf.py` | Validates dependency and library ordering in `config/specs/binaries.conf`. |
| `check-forbidden-flags.py` | Rejects debugging or optimization flags in configuration locations where they are forbidden. Reads `config/specs/testfarm.conf`, the one remaining reason that file is not yet removed -- planned to move to `abibuildbot` (reading `Bconfig.testfarm_flags` there instead) in a future session. |

The commands currently registered in `buildsys/tests.cnf` are the ones invoked by `make tests_buildsys`.
A script merely being present in `abichecks/scripts` does not make it part of that Make target.
It must be registered explicitly in `buildsys/tests.cnf`.

### LibPAW packaging check

`check-libpaw.py` creates a standalone LibPAW archive, extracts it in a temporary location, builds it, and removes the temporary files.
It is registered as check `01` in `libpaw/tests.cnf`.

The `tests_libpaw` target executes it via `run_checks.py`, the same as `tests_abirules`/`tests_buildsys` -- unlike the legacy Perl runner, `run_checks.py` dispatches any series by name and does not hardcode which ones exist.

## Python test-suite registration

The `abirules`, `buildsys`, and `libpaw` directories also contain `__init__.py` files with keywords and `pyscripts` lists.
These files appear to be metadata from an intended or former Python test-suite integration.
No current consumer of these `pyscripts` lists was found elsewhere in the repository.

The Make-driven interface does not read these files.
Its authoritative registration remains `tests.cnf`.
Keep the metadata synchronized when practical, but do not assume that adding a filename to `pyscripts` makes the check executable through either Make or `tests/runtests.py`.

## Adding a new check

To add a check to an existing series:

1. Place the implementation in `abichecks/scripts/`.
2. Make the script return zero on success and nonzero on failure.
3. Keep output concise and include filenames and line numbers when possible.
4. Add a numbered entry to the appropriate `abichecks/SERIES/tests.cnf`.
5. Keep the corresponding `pyscripts` metadata synchronized if the script belongs in that historical list.
6. Run `./config/scripts/makemake` from the repository root.
7. Reconfigure an existing build directory if its generated Makefiles were not refreshed automatically.
8. Run the individual check, the complete series, and `make dist`.

For example, after adding check `22` to the build-system series:

```sh
cd _build/abichecks
make tests_buildsys start=22
make tests_buildsys
cd ..
make dist
```

Choose a stable numeric identifier and preserve the existing grouping convention described by the series README and `tests.cnf` comments.

## Removing or renaming a check

Before removing or renaming a script:

1. Search for the path in `tests.cnf`, `__init__.py`, Make fragments, documentation, and pre-commit configuration.
2. Remove or update the Make-driven registration and any historical `pyscripts` metadata.
3. Remove obsolete input and reference files only after confirming that no other check uses them.
4. Run `makemake` to refresh `config/dist/auto-abichecks.lst` and generated Makefiles.
5. Run the affected series and `make dist` to detect stale distribution entries.

A stale generated distribution list can make `make dist` fail with `No rule to make target ... needed by distdir`.
Always fix the authoritative registration or generator input instead of editing only a generated `Makefile.am`.

## Recommended validation

For changes limited to one maintenance check, use:

```sh
cd _build/abichecks
make tests_SERIES start=NN
make tests_SERIES
```

For changes affecting registration or distributed files, also use:

```sh
cd _build
make dist
```

For release-sensitive build-system or packaging changes, use:

```sh
cd _build
make distcheck
```

Inspect the complete logs even when Make reports success.
Some legacy scripts summarize individual findings in generated report files rather than printing every detail directly to the terminal.
