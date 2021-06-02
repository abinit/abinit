#!/usr/bin/env sh

export ABI_HOME="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo "ABI_HOME set to:\n\t" $ABI_HOME

export PATH=$ABI_HOME/src/98_main/:$PATH     # Do not change this line: path to executable
export ABI_TESTS=$ABI_HOME/tests/            # Do not change this line: path to tests dir
export ABI_PSPDIR=$ABI_TESTS/Psps_for_tests/ # Do not change this line: path to pseudos dir

echo "PATH set to:\n\t" $PATH
echo "ABI_TESTS set to:\n\t" $ABI_TESTS
echo "ABI_PSPDIR set to:\n\t" $ABI_PSPDIR
