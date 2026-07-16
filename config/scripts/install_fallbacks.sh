#!/bin/bash
# Wrapper script placed in the build directory by configure.
# Builds and installs the Abinit fallback libraries.

script_dir=$(dirname "$0")

if [ ! -f "${script_dir}/fallbacks/build-abinit-fallbacks.sh" ]; then
  echo "Error: ${script_dir}/fallbacks/build-abinit-fallbacks.sh not found."
  echo "Make sure you have run configure first."
  exit 1
fi

cd "${script_dir}/fallbacks" && bash build-abinit-fallbacks.sh "$@"
