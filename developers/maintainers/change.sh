#!/usr/bin/env bash
# Copyright (C) 1998-2017 ABINIT group (XG)
# 
# The purpose of this script is to change some
# expression by another in a whole set of files
# change file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file1 
 sed -e 's/temp_/tdepes_/' $file > tmp.file1
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.file1 $file
 echo "file $file written "
done
