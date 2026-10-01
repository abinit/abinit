# ABINIT's CMake build system

This directory contains the modules used by ABINIT's CMake build system.
The top-level `CMakeLists.txt` is the entry point; the files here implement focused parts of configuration, dependency discovery, feature probing, testing, and installation.

## Design philosophy

ABINIT supports both Autotools and CMake.
Autotools remains the build system used to generate the official source distribution, while CMake provides a separate, target-oriented way to configure and build the same source tree.
The two frontends should produce compatible ABINIT executables and configuration symbols, but they do not share all their build metadata.

The CMake implementation follows these principles:

1. Builds must be performed outside the source tree.
2. Compiler and platform capabilities are detected instead of assumed.
3. Each source directory defines a CMake library target and declares its dependencies explicitly.
4. External libraries are represented by imported or namespaced targets such as `BLAS::BLAS`, `MPI::MPI_Fortran`, and `abinit::libxc`.
5. Optional functionality is controlled with cached `ABINIT_*` options.
6. The generated `config.h` preserves the preprocessor interface expected by the ABINIT sources.
7. Installation exports ABINIT targets for downstream CMake projects and also provides `abinit.pc` for `pkg-config` users.

This is a parallel build-system implementation, not a CMake wrapper around `configure` or Make.
Consequently, a change to sources, dependencies, or optional components may have to be represented in both the Autotools metadata and the relevant `CMakeLists.txt` files.

## Configuration flow

Configuration starts in the top-level `CMakeLists.txt` and proceeds approximately as follows:

```text
.current_version
       |
       v
top-level CMakeLists.txt
       |
       +-- cmake/options.cmake
       +-- cmake/compiler_setup.cmake
       +-- dependency discovery and feature probes
       +-- cmake/generate_config_h.cmake --> build/config.h
       |
       +-- shared/CMakeLists.txt --> common support libraries
       +-- src/CMakeLists.txt    --> numbered ABINIT libraries and programs
       |
       +-- cmake/abinit_install.cmake
       +-- cmake/abinit_test.cmake
```

The project version is read directly from `.current_version`.
CMake verifies the C/Fortran interface and generates `fc_mangle.h` with `FortranCInterface`.
It then detects compilers, dependencies, and language features before configuring the top-level `config.h.cmake` template.

The numbered source directories are added in dependency order.
Each directory normally creates a static library, sets a common Fortran module directory, links the targets it depends on, and defines an `abinit::` alias.
The programs in `src/98_main` link to the top of this library graph through `abinit::95_drive`.
That directory also assembles the object libraries into the installable `abinit_lib` archive, whose output name is `libabinit`.

## Building ABINIT with CMake

CMake 3.18 or newer is required.
Always use a separate build directory because `cmake/prevent_build_in_source.cmake` rejects in-source builds.

A typical configuration is:

```sh
cmake -S . -B _build-cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build _build-cmake --parallel
```

If no build type is provided, ABINIT selects `RelWithDebInfo`.
The build writes all static archives under the build-tree `lib/` directory and Fortran module files under `modules/`.
It also enables `compile_commands.json` for editors and language servers.

Use CMake's normal cache-inspection tools to see the complete set of options:

```sh
cmake -S . -B _build-cmake -LAH
cmake --build _build-cmake --target help
```

When changing compilers or making substantial dependency changes, use a fresh build directory rather than reusing an incompatible CMake cache.

## Finding dependencies

The CMake build prefers target-based dependency handling.
Depending on the package, detection uses CMake package configuration files, CMake's built-in `Find*` modules, `pkg-config`, or an ABINIT-specific finder in this directory.

The most useful ways to guide discovery are:

| Mechanism | Use |
|---|---|
| `CMAKE_PREFIX_PATH` | Prefixes containing dependency CMake package files. |
| `PKG_CONFIG_PATH` | Directories containing `.pc` files. |
| `<PackageName>_ROOT` | Root of a particular dependency, as a CMake variable or environment variable. |
| `BLA_VENDOR` | Selection of a particular BLAS/LAPACK implementation. |

The configured compiler environment matters, especially on HPC systems.
Load the desired compiler, MPI, mathematical libraries, NetCDF, and other modules before the first CMake configuration.

Core discovery includes BLAS, LAPACK, MPI, OpenMP, Libxc, HDF5, NetCDF C, and NetCDF Fortran.
Additional modules handle ScaLAPACK, MKL, ELPA, PAPI, XML libraries, Wannier90, Kokkos, YAKL, GPU toolchains, and other optional functionality.
Some external projects can be built as part of the CMake build, but this is opt-in and can require network access or substantially increase build time.

Important options are declared in `cmake/options.cmake`, the top-level `CMakeLists.txt`, and the `build_or_find_*.cmake` modules.
Examples include:

```sh
cmake -S . -B _build-cmake \
  -DABINIT_SCALAPACK_ENABLED=ON \
  -DABINIT_SCALAPACK_FLAVOR=NETLIB \
  -DABINIT_ENABLE_NETCDF_DEFAULT=ON
```

Do not infer the exact option set from this README alone.
The cache output and option declarations are authoritative because supported accelerators and dependencies evolve.

## Feature probes and `config.h`

ABINIT sources use configuration macros regardless of which build frontend is selected.
The CMake path reproduces this contract by collecting results from compiler checks, `try_compile`, and `try_run`, then applying them to the top-level `config.h.cmake` template.

Small probe programs live in `cmake/try_compile/`.
They test Fortran facilities, MPI behavior, OpenMP and offload support, NetCDF interfaces, BLAS behavior, and selected vendor extensions.
Keep probes small and self-contained, and convert their result into the same `HAVE_*` symbol used by the source code.

Cross-compilation requires special care because `try_run` cannot execute target binaries unless an emulator or cached result is supplied.
Prefer a compile-only probe when runtime behavior is not essential.

## Source targets and dependency ordering

The CMake target graph mirrors ABINIT's numbered source hierarchy.
The top-level build adds `shared` first and `src` second.
Their `CMakeLists.txt` files add lower-level directories before higher-level consumers.

A typical source-directory file performs four tasks:

1. Lists the directory's source files in `add_library`.
2. Places generated Fortran modules in `${CMAKE_BINARY_DIR}/modules`.
3. Publishes that module directory to dependent targets.
4. Links lower-level ABINIT and external targets with `target_link_libraries`.

Dependencies should be expressed on targets rather than by adding global include directories or raw linker flags.
Use `PRIVATE`, `PUBLIC`, and `INTERFACE` according to whether a requirement is internal or must propagate to consumers.
Generator expressions should be used when flags apply only to a particular language, configuration, or build/install interface.

## Tests and auxiliary targets

`cmake/abinit_test.cmake` defines Make-like convenience targets that invoke the existing ABINIT test runner.
After building, run for example:

```sh
cmake --build _build-cmake --target check
cmake --build _build-cmake --target test_fast
cmake --build _build-cmake --target tests_in
```

`check` runs the minimal keyword-based test selection from the build tree.
The `test_fast` and `tests_in` targets run the built-in tests used by the other build frontend.
These are custom build targets, not tests registered with CTest, so `ctest` is not the primary interface here.

The top-level `robodoc` target generates the ROBODOC archive inside the build directory and requires the external `robodoc` executable.

## Installation and consumption

Install with:

```sh
cmake --install _build-cmake
```

`cmake/abinit_install.cmake` installs the executables and libraries using `GNUInstallDirs`.
It also generates and installs:

- `abinit-config.cmake` and its version file;
- the exported `abinit::` targets;
- `abinit.pc` for `pkg-config`.

Use `CMAKE_INSTALL_PREFIX` and the standard `CMAKE_INSTALL_BINDIR`, `CMAKE_INSTALL_LIBDIR`, and `CMAKE_INSTALL_INCLUDEDIR` variables to customize the layout.

## Directory map

| Path | Responsibility |
|---|---|
| `options.cmake` | User-facing feature, accelerator, debug, and optimization options. |
| `compiler_setup.cmake` | Compiler-family-specific compile and link flags. |
| `generate_config_h.cmake` | Capability checks and generation of `config.h`. |
| `try_compile/` | Minimal sources used for configure-time feature probes. |
| `build_or_find_*.cmake` | Selection between installed and locally built optional dependencies. |
| `FindWANNIER.cmake` | ABINIT-specific discovery of an installed Wannier90 library. |
| `mkl_setup.cmake`, `scalapack_setup.cmake` | Mathematical-library setup and imported targets. |
| `manage_openmp_target.cmake` | Compiler-specific OpenMP offload flags. |
| `get_mpi_vendor.cmake`, `CheckMPIFeatures.cmake` | MPI implementation and capability detection. |
| `abinit_test.cmake` | Convenience targets for the existing ABINIT test suite. |
| `abinit_install.cmake` | Installation, package export, and `pkg-config` generation. |
| `abinit-config.cmake.in` | Installed CMake package configuration template. |
| `prevent_build_in_source.cmake` | Enforcement of out-of-source builds. |
| `print_target_properties.cmake` | Diagnostic helper for inspecting CMake targets. |

## Maintaining the CMake interface

### Adding or removing a source file

Update the `add_library` or `add_executable` source list in the corresponding source-directory `CMakeLists.txt`.
Make the equivalent change in the Autotools build metadata.
For a new library directory, also add it at the correct dependency position in its parent's `CMakeLists.txt` and link it from its consumers.

Run `./config/scripts/makemake` after changing the source tree or build metadata.
`makemake` does not generate the CMake target graph, but it regenerates the Automake files and distribution lists that ship `CMakeLists.txt` and CMake modules in release archives.

### Adding a dependency or feature

Prefer an imported target supplied by the dependency itself.
If none is available, try `pkg_check_modules(... IMPORTED_TARGET ...)` or add a focused `Find<Package>.cmake` module.
Expose a stable `abinit::<name>` alias when this keeps source-directory logic independent of the discovery method.

Add the user option to `options.cmake` or a focused dependency module.
Perform validation during configuration, set the compatible `HAVE_*` macro, and link the dependency only to targets that need it.
If Autotools supports the same feature, keep option semantics and configuration macros aligned between the two frontends.

### Adding a CMake module or probe

Place reusable modules directly under `cmake/` and probe sources under `cmake/try_compile/`.
Because `cmake` is registered as a detached data directory in `config/specs/buildsys.conf`, run `makemake` and confirm that the new path appears in `config/dist/auto-cmake.lst`.
This step is required for `make dist` to include the file.

## Verification checklist

For changes affecting this build frontend:

1. Configure in a new build directory with the intended compiler and dependency environment.
2. Review the feature summary printed at the end of configuration.
3. Build all default targets.
4. Run `check` and any tests relevant to changed optional features.
5. Run `cmake --install` into a temporary prefix when installation or exports changed.
6. Run `./config/scripts/makemake` and confirm that CMake files are present in the generated distribution lists.
7. Run `make dist` from an Autotools build when distribution content changed.

When modifying shared build behavior, validate both CMake and Autotools rather than assuming success with one frontend proves that the other remains synchronized.
