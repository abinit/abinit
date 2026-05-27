#!/usr/bin/env bash
set -euo pipefail

failed=0

for file in "$@"; do
  if [ ! -f "${file}" ]; then
    continue
  fi

  if grep -Iq . -- "${file}" && grep -n $'\r' -- "${file}" >/dev/null; then
    printf 'CRLF line endings found in %s. Use LF line endings.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
