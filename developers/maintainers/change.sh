#!/usr/bin/env bash
# Copyright (C) 1998-2018 ABINIT group (XG)
# 
# The purpose of this script is to change some
# expression by another in a whole set of files
# change file1 file2 ...

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.file tmp.file2 tmp.file3 tmp.file4 tmp.file5 tmp.file6
 sed -e 's!pp. 477-481.!pp. 477-481. [[cite:Lebedev1999]]!' $file > tmp.file
 sed -e 's!1995, pp. 283-286.!1995, pp. 283-286. [[cite:Lebedev1995]]!' tmp.file > tmp.file2
 sed -e 's!1992, pp. 587-592.!1992, pp. 587-592. [[cite:Lebedev1992]]!' tmp.file2 > tmp.file3
 sed -e 's!1977, pp. 99-107.!1977, pp. 99-107. [[cite:Lebedev1977]]!' tmp.file3 > tmp.file4
 sed -e 's!1976, pp. 10-24.!1976, pp. 10-24. [[cite:Lebedev1976]]!' tmp.file4 > tmp.file5
 sed -e 's!1975, pp. 44-51.!1975, pp. 44-51. [[cite:Lebedev1975]]!' tmp.file5 > tmp.file6
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.file6 $file
 echo "file $file written "
done
