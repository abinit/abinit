#!/bin/bash

# Parse command line options
VERBOSE=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        *)
            echo "Usage: $0 [-v|--verbose]"
            exit 0
            ;;
    esac
done

echo "=== Compiler Information ==="
echo "Fortran compiler: ${FC:-not set}"
echo "C compiler:       ${CC:-not set}"
echo "C++ compiler:     ${CXX:-not set}"

if [ "$VERBOSE" = true ]; then
    echo
    echo "=== Compiler Versions ==="
    if command -v "${FC}" &> /dev/null; then
        echo "Fortran version: $("${FC}" --version 2>&1 | head -1)"
    fi
    if command -v "${CC}" &> /dev/null; then
        echo "C version: $("${CC}" --version 2>&1 | head -1)"
    fi
    if command -v "${CXX}" &> /dev/null; then
        echo "C++ version: $("${CXX}" --version 2>&1 | head -1)"
    fi
    echo
    echo "=== System Information ==="
    echo "Hostname: $(hostname)"
    echo "Architecture: $(uname -m)"
    echo "Operating System: $(uname -s)"
    echo "pkg-config version: $(pkg-config --version 2>&1 || echo 'N/A')"
fi
echo

echo "=== Required Libraries ==="
REQUIRED=("blas" "lapack" "hdf5" "netcdf" "netcdf-fortran")

echo "Checking required library availability using pkg-config..."
for LIB in "${REQUIRED[@]}"; do
    if pkg-config --exists "$LIB"; then
        echo "✅ $LIB"
        echo "   Version: $(pkg-config --modversion "$LIB")"
        echo "   CFLAGS:  $(pkg-config --cflags "$LIB")"
        echo "   LIBS:    $(pkg-config --libs "$LIB")"

        if [ "$VERBOSE" = true ]; then
            echo "   Requires: $(pkg-config --requires "$LIB" 2>/dev/null | tr '\n' ' ' || echo 'none')"
            echo "   Requires-private: $(pkg-config --requires-private "$LIB" 2>/dev/null | tr '\n' ' ' || echo 'none')"
            # Try to get install prefix from pkg-config
            if pkg-config --variable=prefix "$LIB" &> /dev/null; then
                echo "   Prefix: $(pkg-config --variable=prefix "$LIB")"
            fi
        fi
    else
        echo "⚠️  $LIB - NOT FOUND"
    fi
done
echo

# Check optional libraries
echo "=== Optional Libraries ==="
OPTIONAL=("openblas" "mpi" "scalapack" "scalapack-openmpi" "scalapack-mpich" "elpa" "fftw3")

for LIB in "${OPTIONAL[@]}"; do
    if pkg-config --exists "$LIB"; then
        echo "✅ $LIB"
        echo "   Version: $(pkg-config --modversion "$LIB")"
        echo "   CFLAGS:  $(pkg-config --cflags "$LIB")"
        echo "   LIBS:    $(pkg-config --libs "$LIB")"

        if [ "$VERBOSE" = true ]; then
            echo "   Requires: $(pkg-config --requires "$LIB" 2>/dev/null | tr '\n' ' ' || echo 'none')"
            echo "   Requires-private: $(pkg-config --requires-private "$LIB" 2>/dev/null | tr '\n' ' ' || echo 'none')"
            if pkg-config --variable=prefix "$LIB" &> /dev/null; then
                echo "   Prefix: $(pkg-config --variable=prefix "$LIB")"
            fi
        fi
    else
        echo "⚠️  $LIB - not found (optional)"
    fi
done

if [ "$VERBOSE" = true ]; then
    echo
    echo "=== All Available Packages ==="
    pkg-config --list-all 2>/dev/null || echo "Unable to list all packages"
fi
echo

echo "=== Library check complete ==="
exit 0
