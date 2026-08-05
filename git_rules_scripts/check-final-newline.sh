#!/usr/bin/env bash
set -euo pipefail

failed=0

for file in "$@"; do
  if [ ! -f "${file}" ] || [ ! -s "${file}" ]; then
    continue
  fi

  if ! grep -Iq . -- "${file}"; then
    continue
  fi

  last_byte="$(LC_ALL=C tail -c 1 -- "${file}" | od -An -tx1 | tr -d ' \n')"
  if [ "${last_byte}" != "0a" ]; then
    printf '\n' >> "${file}"
    printf 'Added final newline to %s. Stage the fixed file and commit again.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
