#!/bin/bash

#
# Copyright (C) 2010-2018 ABINIT Group (Jean-Michel Beuken)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#set -e

export PATH=/home/buildbot/miniconda3/envs/abimkdocs/bin:/home/buildbot/miniconda3/bin:/home/buildbot/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
#conda info --envs
echo -e "***  Enter in abimkdocs envs...\n"
source activate abimkdocs

echo -e "***  Updating python modules...\n"
pip install -q -r requirements.txt &> /dev/null

echo -e "\n***  Starting SimpleHTTPServer...\n"
cd site
sudo netstat -tapn | grep 8000
#python2 ../web.py &>stderr &
~/bin/web.py &>stderr &
jobs
sleep 1
sudo netstat -tapn | grep 8000

echo -e "\n***  Starting linkchecker...\n"

# BB logfiles={ "link-cleaned" : "linkchecker_ext_wo_false_errors.log" },
# "linkchecker_ext.log" is too big  
echo "cmd : linkchecker -v --no-status --check-extern -o xml --ignore-url=.*fonts.gstatic.com http://localhost:8000/ > ../linkchecker_ext.log"

linkchecker -v --no-status --check-extern --timeout 15 -o xml --ignore-url=.*fonts.gstatic.com http://localhost:8000/ > ../linkchecker_ext.log 2> ../linkchecker_ext.err        

echo -e "\n***  Stopping SimpleHTTPServer...\n"
kill %1
sleep 1
sudo netstat -tapn | grep 8000

cd ..

echo -e "\n***  Removing false errors...\n"

./scripts/bbtools/LinkChecker_rm_false_errors.py linkchecker_ext.log > linkchecker_ext_wo_false_errors.log

rc=`echo $?` 
#echo "exit_status : $rc"
exit $rc
