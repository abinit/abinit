#!/bin/bash

#
# Copyright (C) 2010-2024 ABINIT Group (Jean-Michel Beuken)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#set -e

#conda info --envs

echo -e "***  Remove old mkdocs envs...\n"
conda remove -n mkdocs --all -y
echo -e "***  Create mkdocs envs...\n"
conda create --name mkdocs python=3.7 -y
echo -e "***  Enter in mkdocs envs...\n"
conda activate mkdocs

echo -e "***  Content of requirements.txt\n"
cat requirements.txt

echo -e "***  Updating python modules...\n"
pip install -q -r requirements.txt &> /dev/null

cd abimkdocs_plugin
python setup.py install &> /dev/null
#echo $?
cd ..

echo -e "***  Building the site...\n"

# some patches
sed -i -e "s/google_analytics/#google_analytics/" mkdocs.yml
sed -i -e 's/site_url: https:\/\/docs.abinit.org/site_url: http:\/\/127.0.0.1:8000\//' mkdocs.yml

echo -e "*** Removing modules_with_data...\n"
rm -rf tests/modules_with_data

echo -e "*** Executing mksite.py build...\n"
./mksite.py build
rc=`echo $?`
echo "build_status : $rc"

echo -e "***  Deactivate mkdocs envs...\n"
conda deactivate
conda remove -n mkdocs --all -y

exit $rc
