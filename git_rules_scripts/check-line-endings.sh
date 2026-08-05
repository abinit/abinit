#!/usr/bin/env bash
set -euo pipefail

failed=0

for file in "$@"; do
  if [ ! -f "${file}" ]; then
    continue
  fi

  if grep -Iq . -- "${file}" && grep -n $'\r' -- "${file}" >/dev/null; then
    perl -0pi -e 's/\r\n/\n/g; s/\r/\n/g' -- "${file}"
    printf 'Converted CRLF line endings to LF in %s. Stage the fixed file and commit again.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
