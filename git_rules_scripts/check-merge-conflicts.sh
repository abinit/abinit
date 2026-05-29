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

  if grep -nE '^(<<<<<<<|=======|>>>>>>>)' -- "${file}"; then
    printf 'Possible unresolved merge conflict markers found in %s.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
