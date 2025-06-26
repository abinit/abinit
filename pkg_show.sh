#!/bin/bash

echo "Fortran compiler: ${FC}"
echo "C compiler:       ${CC}"
echo "C++ compiler:     ${CXX}"

# List of libraries to check
LIBRARIES=("netcdf-fortran" "netcdf" "hdf5", "elpa", "fftw3")

echo "Checking library availability and configuration using pkg-config..."

for LIB in "${LIBRARIES[@]}"; do
    if pkg-config --exists "$LIB"; then
        echo "✅ $LIB is available."
        echo "  Version: $(pkg-config --modversion "$LIB")"
        echo "  CFLAGS:  $(pkg-config --cflags "$LIB")"
        echo "  LIBS:    $(pkg-config --libs "$LIB")"
        echo "  PATH:    $(pkg-config --path "$LIB")"
    else
        echo "❌ $LIB is NOT available."
    fi
    echo
done
