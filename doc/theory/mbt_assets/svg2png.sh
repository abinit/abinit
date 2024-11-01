#!/bin/bash
ls *.svg > ./files
sed -i "s/.svg//" ./files
cat ./files | while read filename
do
  convert ./$filename.svg ./$filename.png
done
rm ./files

