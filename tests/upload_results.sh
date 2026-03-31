#!/bin/bash

if [ "$#" != "4" ]; then
 echo " arguments missing : buildername repository sha1 branch"
 exit 99
fi

# check if builder use build dir
Build=0
Path="../"

if test ! -f "${Path}configure.ac"; then
    Path="../../"
    Build=1
    cp testbot.cfg ../../tests/
fi

# ubu_gnu_5.3_openmpi git@gitlab.abinit.org:beuken/abinit8.git ef84dae75e4fb4b9953299f6d3ad9157fff97e86 master

current_builder=$1
current_repository=$2
sha1=$3
current_branch=$4
#current_sha1=`git reflog show --all | cut -d' ' -f1`
current_sha1=`git reflog show --abbrev=8 | cut -d' ' -f1 | head -1`
current_committer=`git config --get remote.origin.url | cut -d: -f2 |cut -d/ -f1`
tag=`git describe --tags --abbrev=0`



ListDirs=`ls -1 -d TestBot_*`


# Send Results to eos
# Since the exit code of the python script that generates the tests summary cannot be read here,
# a ANALYSIS_SUMMARY_FAILED file is created to signal the failure of the analysis. Skip the upload
# if the file is present
if [ ! -f ANALYSIS_SUMMARY_FAILED ]; then
  server="eos"
  remote_base_dir="/data/buildbot_results/bb"
  remote_dir="$remote_base_dir/$current_builder/$current_committer/$current_branch/${current_sha1}"

  ssh $server "test -d $remote_dir"
  rc2=`echo $?`
  if [ "$rc2" ]; then
      ssh $server rm -rf $remote_dir
  else
      remote_dir=${remote_dir}_trial2
  fi

  for i in $ListDirs; do
    ssh $server mkdir -p $remote_dir/$i
    scp $i/results.tar.gz $i/suite_report.html $server:$remote_dir/$i
    ssh $server "cd $remote_dir/$i;tar -mxzf results.tar.gz;chmod -R go+r .;rm -f results.tar.gz"
  done
fi
