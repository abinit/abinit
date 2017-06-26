#!/usr/bin/env bash
# Copyright (C) 1998-2017 ABINIT group (XG)
# 
# The purpose of this script is to change some
# expression by another in a whole set of files
# change file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file
 sed -e 's!varname:!abivarname:!' $file > tmp.file
 sed -e 's!definition:!mnemonics:!' tmp.file > tmp.file2
 echo "changes done "
 # put the modified file at the correct place
 mv $file $file.orig
 mv tmp.file2 $file
 rm tmp.file
 echo "file $file written "
done
