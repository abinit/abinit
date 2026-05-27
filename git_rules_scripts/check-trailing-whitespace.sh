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
    printf 'Trailing whitespace found in %s. Remove spaces or tabs at line ends.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
