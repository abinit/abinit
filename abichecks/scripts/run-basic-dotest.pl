# This file will run one abinit built-in test in a separate process

# Copyright (C) 1999-2020 ABINIT group (LSi)
# This file is distributed under the terms of the
# GNU General Public License, see ~ABINIT/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~ABINIT/Infos/contributors .
#
# Usage :
# This script is intended to be called from script run-basic-tests[.pl]
# The latter starts the test in background mode and waits for
# completion. When done, this script creates a test#.end file if required
#
# Arguments :
#   $TestN = test number (fast, v1, v5, bigdft, etsf_io, libxc, wannier90)
#   $ENDtest = kind of criterion for test completion; supported values are:
#	ps = unix ps command
#	pstat = Windows NT ResKit pstat command
#	numeric = test#.end file with loop limit
#   $SLASH = subdirectory delimitor in file paths
#   $SUFXstyle = defines suffixes style for binaries, ...;
#		supported values are: unix, DOS, vms
#
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator
$debug = 0;             # debug level
#
$timeout = $ENV{'timeout'};
$timeout = '0' if ($timeout eq '');
#
($TestN,$ENDtest,$SLASH,$SUFXstyle) = @ARGV;
if ($SUFXstyle eq 'unix') {
# unix suffixes:
	$XSUFX = '';		# no special suffix for binary module
	}
elsif ($SUFXstyle eq 'DOS') {
# DOS suffixes:
	$XSUFX = '.exe';	# suffix for binary module
	}
elsif ($SUFXstyle eq 'vms') {
# DOS suffixes:
	$XSUFX = '.exe';	# suffix for binary module
	}
# make sure "files" file exist
if (! -r "test$TestN.files") {
	print "Test number missing or invalid: $TestN";
	exit 8;
	}
print "run-basic-dotest.pl: Starting built-in test $TestN" if ($debug >= 1);
$abinit = &transpath("../../../src/98_main/abinit$XSUFX");
$timeout_cmd = &transpath("../../Timeout/timeout$XSUFX");
if ($SUFXstyle eq 'vms') {
  $cmd = "pipe run [-]abinit < test$TestN.files > test$TestN.log";
  }
else {
  $cmd = "$abinit < test$TestN.files > test$TestN.log";
  }
if ($timeout ne '0') {
  $cmd = "$timeout_cmd $timeout $cmd";
  }
print $cmd if ($debug >= 1);
$rc = system($cmd);		# start abinit
$_ = $ENDtest;
# check numeric value for test completion
if (! /[^0-9]/) {
	open (END,">test$TestN.end");	# create a file to signal end of test to caller
	close (END);
	}
exit;
# ****************************************
sub transpath {
	local($path) = @_;
#
# purpose: translate unix-like path of file to DOS-like according to $SLASH
# argument:
#	$path = path to be translated
# output: this subroutine acts as a function and return the path
# according to host conventions
	
	$path =~ tr/\057/\\/ if ($SLASH eq '\\');
	return $path;
	}
