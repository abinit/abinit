#!/bin/bash

for f in `find src -name '*.F90'`; do
  src_string=`echo ${f#src/} | sed -e 's:/:,:g'`
  pcregrep '^!!' $f | \
    pcregrep -M '(!!\*\*\*\*f|PARENTS\n!!\n!! CHILDREN)' | \
    pcregrep -M '!!\*\*\*\*f\* [A-Za-z0-9_/-]*\n!! PARENTS' | \
    grep -v PARENTS | \
    sed -e 's,.*/,,' | \
    awk "{printf \"%s,%s\n\",\"$src_string\",\$1}"
done | sort | tee abinit-orphans-7.9.1.log
