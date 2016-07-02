# extract differences with context from fldiff output
file='fldiff.report'
if test X$1 != X; then
file=$1
fi
if  test ! -e $file; then
echo File $file does not exist
exit 16
fi
echo Case line# context
awk 'BEGIN {num="??"} /^ Case/ {num=substr($1,6)} /^[0-9]/ {if ($2 != "" && $2 != "Delivered" && $4 != "WARNINGs") print num,$0}' $file | sort +2.0 +3.0
