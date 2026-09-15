#!/bin/sh

# Build the ABINIT source-code reference in a reproducible way.
# Set DOXYGEN_WARNINGS_AS_ERRORS=1 in CI to reject builds with Doxygen warnings.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$script_dir"

if ! command -v doxygen >/dev/null 2>&1; then
    echo "[mkdoxygen] Error: doxygen is not available in PATH." >&2
    echo "[mkdoxygen] Install Doxygen and Graphviz, then run ./mkdoxygen.sh again." >&2
    exit 127
fi

if ! command -v dot >/dev/null 2>&1; then
    echo "[mkdoxygen] Error: Graphviz 'dot' is not available in PATH." >&2
    echo "[mkdoxygen] Doxyfile enables dependency graphs, so Graphviz is required." >&2
    exit 127
fi

if [ ! -f Doxyfile ] || [ ! -f .current_version ]; then
    echo "[mkdoxygen] Error: Doxyfile or .current_version is missing from the repository root." >&2
    exit 2
fi

ABINIT_DOXYGEN_VERSION=$(tr -d '[:space:]' < .current_version)
if [ -z "$ABINIT_DOXYGEN_VERSION" ]; then
    echo "[mkdoxygen] Error: .current_version does not contain an ABINIT version." >&2
    exit 2
fi
export ABINIT_DOXYGEN_VERSION

echo "[mkdoxygen] Building ABINIT $ABINIT_DOXYGEN_VERSION source-code reference."
echo "[mkdoxygen] Removing previous generated documentation files."
rm -rf -- doxygen_docs
rm -f -- doxygen.err abinit.tag

if doxygen -q Doxyfile; then
    doxygen_status=0
else
    doxygen_status=$?
fi

warning_count=0
if [ -f doxygen.err ]; then
    warning_count=$(awk 'NF { count++ } END { print count + 0 }' doxygen.err)
fi

echo "[mkdoxygen] Doxygen exit status: $doxygen_status"
echo "[mkdoxygen] Diagnostic lines: $warning_count (see doxygen.err)"

if [ "$doxygen_status" -ne 0 ]; then
    exit "$doxygen_status"
fi

if [ ! -f doxygen_docs/html/index.html ] || [ ! -f doxygen_docs/xml/index.xml ] || [ ! -f abinit.tag ]; then
    echo "[mkdoxygen] Error: Doxygen did not produce all expected HTML, XML, and tag outputs." >&2
    exit 1
fi

if [ "${DOXYGEN_WARNINGS_AS_ERRORS:-0}" = "1" ] && [ "$warning_count" -ne 0 ]; then
    echo "[mkdoxygen] Error: warnings are treated as errors by DOXYGEN_WARNINGS_AS_ERRORS=1." >&2
    exit 1
fi

echo "[mkdoxygen] Documentation available at doxygen_docs/html/index.html"
echo "[mkdoxygen] Machine-readable output available at doxygen_docs/xml/index.xml and abinit.tag"
