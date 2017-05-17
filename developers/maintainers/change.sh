#!/usr/bin/env bash
# Copyright (C) 1998-2017 ABINIT group (XG)
# 
# The purpose of this script is to change some
# expression by another in a whole set of files
# change file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file1 tmp.file2 tmp.file3
 sed -e 's/#%%<BEGIN TEST_INFO>/AAAA/' $file > tmp.file1
 sed -e '/AAAA/a\\#%%<BEGIN TEST_INFO>' tmp.file1 > tmp.file2
 sed -e 's!AAAA!## After modifying the following section, one might need to regenerate the pickle database with runtests.py -r!' tmp.file2 > tmp.file3
#sed -e 's!#%%<BEGIN TEST_INFO>!## After modifying the following section, one might need to regenerate the pickle database with runtests.py -r \\#%%<BEGIN TEST_INFO>!' $file > tmp.file
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.file3 $file
 echo "file $file written "
done
