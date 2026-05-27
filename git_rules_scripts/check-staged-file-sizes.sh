#!/usr/bin/env bash
set -euo pipefail

max_file_size_bytes="${MAX_FILE_SIZE_BYTES:-11534336}"
oversized=0

while IFS= read -r -d '' file; do
  if ! git cat-file -e ":${file}" 2>/dev/null; then
    continue
  fi

  size="$(git cat-file -s ":${file}")"
  if [ "${size}" -gt "${max_file_size_bytes}" ]; then
    printf 'File exceeds 11 MiB limit: %s (%s bytes)\n' "${file}" "${size}"
    oversized=1
  fi
done < <(git diff --cached --name-only -z --diff-filter=ACMRT)

if [ "${oversized}" -ne 0 ]; then
  printf 'Commit rejected: remove or reduce staged files larger than 11 MiB.\n'
  exit 1
fi

printf 'All staged files are within the 11 MiB limit.\n'
