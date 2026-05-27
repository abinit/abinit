#!/usr/bin/env bash
set -euo pipefail

failed=0

for file in "$@"; do
  if [ ! -f "${file}" ] || [ ! -x "${file}" ]; then
    continue
  fi

  first_two_bytes="$(LC_ALL=C head -c 2 -- "${file}" | od -An -tx1 | tr -d ' \n')"
  if [ "${first_two_bytes}" != "2321" ]; then
    printf 'Executable file %s is missing a shebang line. Add one, for example #!/usr/bin/env bash.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
