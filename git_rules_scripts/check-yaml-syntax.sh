#!/usr/bin/env bash
set -euo pipefail

failed=0

for file in "$@"; do
  case "${file}" in
    *.yaml|*.yml) ;;
    *) continue ;;
  esac

  if [ ! -f "${file}" ]; then
    continue
  fi

  if command -v ruby >/dev/null 2>&1; then
    if ! ruby -e 'require "yaml"; YAML.load_file(ARGV.fetch(0))' "${file}"; then
      printf 'Invalid YAML syntax in %s.\n' "${file}"
      failed=1
    fi
  elif command -v python3 >/dev/null 2>&1; then
    if ! python3 -c 'import sys, yaml; yaml.safe_load(open(sys.argv[1], encoding="utf-8"))' "${file}"; then
      printf 'Invalid YAML syntax in %s.\n' "${file}"
      failed=1
    fi
  else
    printf 'Cannot check YAML syntax for %s: install ruby or python3 with PyYAML.\n' "${file}"
    failed=1
  fi
done

exit "${failed}"
