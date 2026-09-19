# ABINIT build-system configuration

This directory contains the inputs and generators used to maintain ABINIT's Autotools build system.

The central principle is that maintainers describe the source tree and build rules in `config/`, then run `config/scripts/makemake` to regenerate the Autotools files.
Files such as the top-level `Makefile.am`, `configure.ac`, and the `auto-*` files are generated outputs or contain generated sections.
Their headers indicate whether they may be edited directly.
Do not fix a generated file without also changing the configuration or generator that produces it, because the next `makemake` run will overwrite the fix.

## Directory map

The most important parts of `config/` are:

| Path | Purpose |
|---|---|
| `specs/` | Declarative descriptions of the source tree, options, dependencies, and files to ignore. |
| `scripts/` | Programs that generate Makefiles, M4 macros, source lists, and other build-system files. |
| `makefiles/` | Hand-maintained fragments inserted into generated `Makefile.am` files. |
| `m4/` | Hand-maintained and generated Autoconf macros. Generated macros use the `auto-` prefix. |
| `dist/` | Generated lists of files that must be included in source distributions. |
| `gnu/` | Auxiliary files installed by Autoconf and Automake. |
| `templates/` | Templates consumed by build-system generators. |
| `examples/`, `hints/`, `optim/` | Example configurations and compiler/platform settings. |

The top-level `scripts/` directory is different from `config/scripts/`.
The former contains utilities distributed to ABINIT users, while the latter implements the build-system machinery.

## Why `makemake` is important

Run `makemake` from the top-level ABINIT source directory:

```sh
./config/scripts/makemake
```

`makemake` synchronizes several representations of the project.
Among other operations, it:

1. Reads the files under `config/specs/`.
2. Regenerates source and dependency information.
3. Regenerates `config/dist/auto-*.lst` distribution lists.
4. Regenerates the top-level and recursive `Makefile.am` files.
5. Regenerates `config/m4/auto-*.m4` and related Autoconf input.
6. Runs the Autotools programs that produce `configure` and `Makefile.in` files.

An executable Git pre-commit hook is normally required.
Install the hooks from the top-level directory, where `.pre-commit-config.yaml` is located:

```sh
python3 -m pip install pre-commit
pre-commit install
```

For exceptional or automated maintenance where hooks cannot be installed, use:

```sh
./config/scripts/makemake --no-precommit-hooks
```

Use `./config/scripts/makemake --help` to see partial-generation options.
A full run is preferable before committing build-system changes because partial runs can leave generated files out of sync.

Review all generated changes after running the command.
Do not assume that every changed generated file is unrelated merely because it was not edited by hand.

## How directories are registered

`config/specs/buildsys.conf` divides the repository into logical blocks.
Its main block types are:

- `master`: the core tree controlled by the top-level build system.
- `subsystem`: a component with its own build system that may be attached to or detached from the core.
- `data`: passive content that is distributed but not recursively built.

For example, the documentation block registers several data directories:

```ini
[doc]
type = data
mode = detached
subdirs = doc abimkdocs abinit_theme scripts
```

The order of directories in `SUBDIRS` is generated from these blocks and their priorities.
Only directories participating in the recursive build appear in `SUBDIRS`.
Passive data directories can be included in a source archive without being configured or built.

`config/specs/junk.conf` defines files and directories that the file-list generator must ignore.
This prevents build products, caches, editor files, and other local artifacts from entering generated distribution lists.

## File lists and `make dist`

The distribution path is:

```text
config/specs/buildsys.conf + config/specs/junk.conf
                         |
                         v
          config/scripts/make-file-lists
                         |
                         v
              config/dist/auto-*.lst
                         |
                         v
          config/scripts/make-makefiles-top
                         |
                         v
                    Makefile.am
                         |
                         v
                     make dist
```

`make-file-lists` walks the registered directories, applies the exclusions from `junk.conf`, and writes one `config/dist/auto-DIRECTORY.lst` file per directory.
`make-makefiles-top` incorporates these lists into the generated top-level `Makefile.am` as `EXTRA_DIST` entries.

Automake's `make dist` then:

1. Creates a temporary `abinit-VERSION/` distribution directory.
2. Recursively collects files managed by Automake subdirectories.
3. Copies files named by `EXTRA_DIST`, including the generated `auto-*.lst` content.
4. Runs the project `dist-hook`.
5. Creates the source archive, normally `abinit-VERSION.tar.gz`.
6. Removes the temporary distribution directory.

The `dist-hook` writes `.tarball-version` into the archive and checks that the configured version agrees with `config/scripts/git-version-gen`.
This protects releases from being produced with stale version metadata.

`make distcheck` goes further.
It creates the archive, unpacks it into a separate tree, configures and builds it, runs the test targets, checks installation and uninstallation, and verifies that the distributed source is self-contained.
Use it for release validation whenever the required dependencies and execution time are available.

## Adding a file

### Adding a file inside a registered directory

For a normal source or data file under a directory already listed in `buildsys.conf`:

1. Add the file to the repository.
2. Check that it is not accidentally matched by `config/specs/junk.conf`.
3. Run `./config/scripts/makemake`.
4. Confirm that the appropriate `config/dist/auto-*.lst` contains the file.
5. Review the regenerated `Makefile.am` and related files.
6. Run the relevant build or tests, followed by `make dist` when distribution content changed.

Do not add the file directly to an `auto-*.lst` file.
Those lists are regenerated.

### Adding a new top-level data directory

If the directory contains passive files that should be shipped but not built:

1. Add it to the appropriate `subdirs` entry in `config/specs/buildsys.conf`, or create a well-described data block.
2. Add directory-specific exclusions to `config/specs/junk.conf` if necessary.
3. Run `makemake`.
4. Check that `config/dist/auto-NEWDIR.lst` was generated and is referenced by `Makefile.am`.
5. Run `make dist` and inspect the archive.

If the new directory contains a component that must be configured and built, it may require a `master` or `subsystem` entry, a priority, Makefile-generation support, and possibly Autoconf macros.
Use an existing component with the same role as the model rather than registering it as passive data.

### Adding exceptional top-level files

Some files do not naturally belong to a registered directory.
Hand-maintained additions to the top-level distribution are declared in `config/makefiles/top.am`.
The contents of this fragment are copied into the generated top-level `Makefile.am`.

Add an exceptional file to `config/makefiles/top.am`, run `makemake`, and verify the resulting `EXTRA_DIST` entry.
Avoid listing an entire directory there if it should instead be described in `buildsys.conf` and tracked through an automatically generated file list.

## Removing a file or directory

### Removing files from an existing registered directory

1. Remove the files from the repository.
2. Search for code, documentation, tests, packaging recipes, and configuration that refer to them.
3. Run `./config/scripts/makemake` immediately.
4. Confirm that the removed paths disappeared from the corresponding `config/dist/auto-*.lst`.
5. Run `make dist` to detect stale distribution references.

If the generated list is not refreshed, `make dist` may fail when Automake tries to copy a file that no longer exists.

### Removing a subdirectory while keeping its registered parent

A subdirectory such as `scripts/example-tools` does not normally need a change in `buildsys.conf` when its parent `scripts` remains registered.
Remove the subdirectory, update references, run `makemake`, and confirm that its entries disappeared from `config/dist/auto-scripts.lst`.

Also inspect custom packaging targets in `config/makefiles/top.am`.
Targets such as `dist-lite`, `binary_package`, and `bin_prep` may manipulate particular directories independently of the normal Automake file lists.

### Removing an entire registered directory

1. Remove the directory from `config/specs/buildsys.conf`.
2. Remove directory-specific rules from `config/specs/junk.conf` when they are no longer useful.
3. Search `config/makefiles/`, `configure.ac`, CMake files, documentation, tests, and scripts for explicit references.
4. Remove custom `EXTRA_DIST`, install, cleanup, archive, or `mv` commands that expect the directory.
5. Delete the directory.
6. Run `makemake` and verify that its `config/dist/auto-DIRECTORY.lst` and generated Makefile references disappear.
7. Run the normal tests and `make dist`.
8. Run `make distcheck` for release-sensitive changes.

Removing a directory from `SUBDIRS` or `EXTRA_DIST` alone is not sufficient when custom packaging targets still name it.

## Verifying distribution changes

After adding or removing distributed content, a useful verification sequence is:

```sh
./config/scripts/makemake
make dist
tar -tzf abinit-*.tar.gz | less
```

Check that newly added files are present, removed files are absent, and no generated or local files leaked into the archive.
Use an explicit archive name instead of the wildcard if several ABINIT archives are present.

For stronger validation, run:

```sh
make distcheck
```

Always inspect `git diff` and `git status` before committing.
Generated distribution lists, Makefile inputs, and Autotools outputs are part of the same logical change and should remain synchronized.

## Troubleshooting

### `make dist` reports a missing file

The most common cause is a stale `config/dist/auto-*.lst` or a stale explicit entry in `config/makefiles/top.am`.
Search both the generated and source configuration:

```sh
rg 'missing/path' config Makefile.am
```

Correct the authoritative configuration and run `makemake` again.

### A file is missing from the archive

Check that its top-level directory is registered in `config/specs/buildsys.conf` and that `config/specs/junk.conf` does not exclude it.
Then inspect the corresponding `config/dist/auto-*.lst` and the generated `Makefile.am`.

### A removed directory reappears in generated files

Search the generator inputs rather than editing the generated output:

```sh
rg 'directory-name' config/specs config/makefiles config/scripts configure.ac
```

Explicit rules in `config/makefiles/top.am` are especially important for legacy binary-package and reduced-distribution targets.
