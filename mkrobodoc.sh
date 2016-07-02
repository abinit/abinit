#!/bin/sh
# This script generates the abinit documentation in ROBODOC format in the directory tmp-robodoc
echo "Will generate ROBODOC documentation in tmp-robodoc (requires robodoc)"

rm -rf tmp-robodoc robodoc-html && mkdir tmp-robodoc
cp -rf ./src/[0-9]* tmp-robodoc
cp ./config/robodoc/robodoc-html.rc tmp-robodoc/robodoc.rc
cd tmp-robodoc && rm */*.in && rm */interfaces* && robodoc > ../robodoc.log 2> ../robodoc.err
exit_status=`cat ../robodoc.err | wc -l`
if test $exit_status -ne 0 ; then 
  cat ../doc/developers/robodoc.doc.txt >> robodoc.err
  cat ../robodoc.err
fi

echo "Exit status: " $exit_status
exit $exit_status
