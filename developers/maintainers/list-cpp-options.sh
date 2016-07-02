#!/bin/sh

cpp_options=`grep -rh '^#' src/*/*.F90 | sed -e 's/#[ ]*//;s/#//' | grep '^if' | sed -e 's/if[^ ]* //;s/[!()&|]//g;s/defined / /g;' | perl -e 'while ( <> ){ s/[\s]+/ /g; s/ /\n/g; print $_; }' | sort -u`

cat <<EOF
Preprocessing options used in ABINIT
------------------------------------

${cpp_options}

EOF

for opt in ${cpp_options}; do
   opt_files=`grep -lr "${opt}" src/*/*.F90`
   opt_count=`echo "${opt_files}" | wc -l`
   opt_uline=`echo "${opt}" | sed -e 's/./-/g'`

   echo "${opt}"
   echo "${opt_uline}"
   echo ""
   echo "Used ${opt_count} times"
   echo ""
   if test "${opt}" != "HAVE_CONFIG_H"; then
      echo "${opt_files}"
      echo ""
   fi
done

cat <<EOF
Miscellaneous
-------------

Files with Fortran comments on lines containing CPP directives:

EOF

grep -lr '^#else.*!' src/*/*.F90 | sed -e 's/^/ * /'
grep -lr '^#endif.*!' src/*/*.F90 | sed -e 's/^/ * /'
echo ""
