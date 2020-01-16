#!/usr/bin/env bash
# Copyright (C) 1998-2019 ABINIT group (XG)
# 
# The purpose of this script is to suppress lines that contain some
# expression in a whole set of files
# suppress file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file 
 sed -n '/src2tex/!p' $file > tmp.file
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.file $file
 echo "file $file written "
done
