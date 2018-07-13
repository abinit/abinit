#!/usr/bin/env bash
# Copyright (C) 1998-2018 ABINIT group (XG)
# 
# The purpose of this script is to change some
# expression by another in a whole set of files
# change file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file tmp.file2 tmp.file3 tmp.file4
 sed -e 's!~abinit/doc/users/acknowledgments.html!https://docs.abinit.org/theory/acknowledgments!' $file > tmp.file
 sed -e 's!https://www.abinit.org/about/?text=acknowledgments!https://docs.abinit.org/theory/acknowledgments!' tmp.file > tmp.file2
 sed -e 's!http://www.abinit.org/about/?text=acknowledgments!https://docs.abinit.org/theory/acknowledgments!' tmp.file2 > tmp.file3
 sed -e 's!~abinit/doc/users/acknowledgments.htm!https://docs.abinit.org/theory/acknowledgments!' tmp_file3 > tmp.file4
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.file4 $file
 echo "file $file written "
done
