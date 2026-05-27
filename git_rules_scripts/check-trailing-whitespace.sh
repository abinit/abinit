#!/usr/bin/env bash
set -euo pipefail

failed=0

for file in "$@"; do
  if [ ! -f "${file}" ]; then
    continue
  fi

  if ! grep -Iq . -- "${file}"; then
    continue
  fi

  if grep -n '[[:blank:]]$' -- "${file}"; then
    perl -0pi -e 's/[ \t]+$//mg' -- "${file}"
    printf 'Removed trailing whitespace from %s. Stage the fixed file and commit again.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
