# ABINIT source-code reference

This site is the API and source-code reference generated from ABINIT's Fortran and C sources.

## Browse the code

- Use **Namespaces** to find Fortran modules and their procedures.
- Use **Files** to inspect source files and their documented entities.
- Use **Directories** to follow the organization of the source tree.
- Use the search box to locate a symbol by name.

The reference covers the main sources in `src`, the common support code in `shared/common/src`, and LibPAW in
`shared/libpaw/src`.

## Related documentation

The [ABINIT documentation](https://docs.abinit.org/) contains the user guide, tutorials, input-variable reference, and
developer documentation.

## Build locally

From the repository root, run `./mkdoxygen.sh`.
The generated entry point is `doxygen_docs/html/index.html`, and diagnostics are written to `doxygen.err`.
Structured output for documentation tools and AI agents is available in `doxygen_docs/xml`, with cross-reference data
in `abinit.tag`.
