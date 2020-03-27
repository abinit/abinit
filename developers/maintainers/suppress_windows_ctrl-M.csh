#!/bin/csh
#
# Copyright (C) 2001-2020 ABINIT group (Lsi,XG)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#
# For a DOS-formatted file, each line end with a
# CR-NL (Carriage Return-New Line), while with
# Unix-formatted files, only the NL is used. 
# Under Unix, some commands recognize the CR-NL as being
# simply a new line (e.g. vi), while some other do not (e.g. diff).
# This script allows to suppress the CR at
# the end of each line of a DOS-formatted file, so as to
# make it a Unix-formatted file.

# Usage :
# suppress_windows_ctrl-M.csh file1 file2 ...

# Count the number of files
set number_arg=$#argv
set index=1
echo "the following files will be treated :"
while ($index <= $number_arg)
  echo " $argv[$index] "
  @ index++
end

echo " number of files = $number_arg "

# Set up a loop on all the files
set index=1
while ($index <= $number_arg)

  # Examine each file
  set file=$argv[$index]
  echo " file $file is now treated "
  rm -f tmp.file

  # Here is the suppression of the CR
  awk '/\r$/ {sub("\r$","");print $0}' $file > tmp.file 
  echo " changes done "
  mv $file $file.OLD
  mv tmp.file $file
  echo " file written "

  # Update index
  @ index++
  echo " index is updated to $index "

# End of the loop on the files
end


