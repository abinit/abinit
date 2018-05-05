#!/usr/bin/env bash
# Copyright (C) 1998-2018 ABINIT group (XG)
# 
# The purpose of this script is to change some
# expression by another in a whole set of files
# change file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file  tmp_file2
 sed -e 's!53_spacepar!54_spacepar!' $file > tmp.file
 sed -e 's!54_abiutil!55_abiutil!' tmp_file > tmp.file2
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.file2 $file
 echo "file $file written "
done
