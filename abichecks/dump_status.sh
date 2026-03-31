#!/bin/bash

for i in $(ls -1 *.status); do
  label=`echo $i | cut -d. -f1`
  if [ "$label" == "luo" ]; then label='Less Used Options'; fi
  echo ${label^^}
  echo "------------------"
  cat $i
  echo
done

#exit 0

echo "ABIRULES"
echo "------------------"
cat abirules/report.log
echo

echo "BUILDSYS"
echo "------------------"
cat buildsys/report.log
echo


# generate abichecks/testbot_summary.json
rc=`egrep -i FAIL *.status */report.log`

if [ "rc" == "" ]; then
   echo "{ \"debug\": \"succeeded\" }" > testbot_summary.json
else
   echo "{ \"debug\" : \"failed\" }" > testbot_summary.json
fi
