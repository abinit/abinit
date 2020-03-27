#!/usr/bin/perl
#
# Script for running various tests of abinit, mrgddb, anaddb, etc ...
# It is designed to test the features introduced in the development
# of the successive versions of the ABINIT code.
# Several suites of tests are defined by the tests.cnf file in the fast,
# v1, v2, ... directories.
# Before running this script, read README in this directory.

# Copyright (C) 2004-2020 ABINIT group (LSi)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .

# Usage :
# unix C-shell: run-standard-tests machine_name top_bindir top_testdir [serie] [first_test# [last_test# | end] ] >& log_file
# unix bash: run-standard-tests machine_name top_bindir top_testdir [serie] [first_test# [last_test# | end] ] >log_file 2>&1
# Windows DOS box: perl run-standard-tests.pl machine_name top_bindir top_testdir [serie] [first_test# [last_test# | end] ] >log_file
# MacOS (likely...) : perl -f run-standard-tests.pl machine_name top_bindir top_testdir [serie] [first_test# [last_test# | end] ] >log_file
# 	where serie can be fast, v1, v2, v3, v4, v5, v6, v67mbpt, v7, physics, cpu or tutorial; if omitted, the tests
# configuration file "tests.cnf" is read from the current directory.
# The top_testdir must be an absolute path.
#
# For example, under unix C-shell, on the machine fock :
#  run-standard-tests fock ../src/98_main /home/me/abinit/test fast  >& log_file        run all tests in fast on hilbert
#  run-standard-tests hilbert ../src/98_main /home/me/abinit/test v1 1  >& log_file        run only test 1 in v1 on hilbert
#  run-standard-tests hilbert ../src/98_main /home/me/abinit/test v2 3 end  >& log_file    run tests from case_3 to the end in v2 on hilbert
#  run-standard-tests hilbert ../src/98_main /home/me/abinit/test tutorial 8 22  >& log_file    run tutorial tests between 8 and 22 included.
#

# Perl path should be defined either by the environmental variable PERL in DOS
# format or in the usual PATH. This is required for cygwin or DOS.
$PERL = $ENV{'PERL'};
$PERL = 'perl' if ($PERL eq '');

# Timeout support
# This is required for nightly builds, because e.g. binaries built with ifort 10.0 make
# the process hang forever when they crash.
$timeout = $ENV{'timeout'};
$timeout = '0' if ($timeout eq '');

$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator
@Series = ('abirules','buildsys');
$TestDir = '.';		# default is current directory
# abinit processing options:
$optKeepfn4 = 1;	# keep files fn4* after test completion
$optistwfk1 = 2;	# for cpu tests, set istwfk=1
$optistwfk11 = 3;	# for cpu tests, set istwfk=1 1
# diff/fldiff processing options:
$fldNorun = 0;	# don't run the fldiff script
$fldRun = 1;	# run the fldiff script (default)
$fldDOS = 2;	# run diff/fldiff on the file tNNNo_DOS
$fldxml = 3;	# run diff/fldiff also on the xml file tNNN.o.cml
$fldDS2 = 4;	# run diff/fldiff on the file tNNNo_DS2_DOS
$fldlog = 5;	# log to .out file and run diff/fldiff on it
$fldlgx = 6;	# run diff/fldiff on the xml file tNNN.o_LOG.xml
$fldMDF = 7;	# run diff/fldiff on the MDF files produced by the BSE code
$deffldopt = '-ignore -ignoreP';	# default options for fldiff script
# ??
$debug = 0;

open(STDOUT,">-"); select(STDOUT); $| = 1;      # make unbuffered

$abinit_srcdir = $ENV{"abinit_srcdir"};
$abinit_builddir = $ENV{"abinit_builddir"};
$abinit_nightlydir = $ENV{"nightly_bindir"};
$abinit_bindir = $ENV{"abinit_bindir"};
$abinit_inpdir = $ENV{"abinit_inpdir"};
$abinit_outdir = $ENV{"abinit_outdir"};
$top_testdir = $ENV{"abinit_inpdir"};
$java_runner = $ENV{"run_java"};

if ($ARGV[0] eq '') {
  print 'an argument must be provided giving machine name';
  exit 16;
  }
$MACH = $ARGV[0];
print "Machine $MACH"  ;
#$abinit_bindir = $ARGV[1];
print "Bindir $abinit_bindir"  ;
#$top_testdir = $ARGV[2];
print "Top_testdir $top_testdir"  ;
$argptr = 1;
# check first argument against some tests serie
foreach $serie (@Series) {
	if ($ARGV[$argptr] eq $serie) {
		$TestDir = $ARGV[$argptr];
		$argptr ++;
		last;
		}
	}
print "Testdir  $TestDir"  ;
# set default first and last test number
$FIRST = 'first';
$LAST = 'end';
# check for first_test# specified
if ($ARGV[$argptr] ne '') {
	$FIRST = $ARGV[$argptr];
# if last_test# is omitted, run only first
	if ($ARGV[$argptr+1] eq '') {
		$LAST = $FIRST;
		}
# check for last_test specified
	elsif ($ARGV[$argptr+2] eq '') {
# change last test number if not end
		$LAST = $ARGV[$argptr+1] if ($ARGV[$argptr+1] ne 'end');
		}
	else {
		print 'Too many arguments';
		exit 12;
		}
	}
# check operating system and set environment accordingly
$UNIXSLASH = '/';	# unix subdirectory delimitor in file paths (octal 057)
$DOSSLASH = '\\';	# subdirectory delimitor in DOS file paths (\ escaped)
$MACSLASH = ':';  # subdirectory delimitor in MAC file paths
$MULTI = '*';
# fetch some environment variables and try unix 'uname' command (Bourne shell)
$OSname = $ENV{'OS'};		# OS name on Windows systems
$OStype = $ENV{'OSTYPE'};	# defined on some systems as OPENVMS, mingw/msys
# Use internal variable of Perl $^O
if ($^O eq 'MacOS') {
  $OStype = 'MacOS';
  $MULTI = "\xc5";
  }
elsif ($OStype ne 'OPENVMS') {
  $unamerc = open(IP,"uname -a |");
  $_ = <IP>;
  chop $_;
  ($UNkernel,$wd2,$wd3,$wd4,$wd5,$wd6,$wd7,$wd8,$wd9) = split(' ',$_);
  close IP;
  }
print "uname: =$OSname=$OStype=$UNkernel=$unamerc= " if ($debug >= 1);
# Check for Windows 2000/XP (DOS box, cygwin, mingw, PGI Workstation, ...)
# cygwin, mingw, PGI Workstation are unix bash shells. However, older versions
# of the latter behave as a DOS one on some aspects
if ($OSname eq 'Windows_NT') {	# check Windows environmental variable
# When 'uname' command is unknown, the return code is undetermined.
# drop version code from kernel name
# by JMB : add support of mixed env cygwin + mingw compilation : MINGW64_NT
  $UNkernel = 'MINGW64_NT' if ( $UNkernel eq 'CYGWIN_NT-6.1-WOW64');
  $UNkernel = 'MINGW64_NT' if ( $UNkernel eq 'CYGWIN_NT-5.2-WOW64');
  $UNkernel =~ m/(.*)(-)([0-9]+\.[0-9]+)/;
  $UNkernel = $1 if ($& ne '');
  if ($OStype eq '' && $UNkernel eq '') {
    $OStype = 'DOS';        # set OStype
    $SLASH = $DOSSLASH;		# subdirectory delimitor in DOS file paths
    $COPY_COMMAND = 'copy';		# DOS copy file command
    $DIFF_COMMAND = 'fc /w /l /n';	# DOS difference command
# since unlink <file*> fails on some versions of perl for DOS, let's use del
    $ERASE_COMMAND = 'del /q';	# DOS delete file command
    $ERASE_RECURSIVE = 'del /s /q';	# NT DOS command to recursively erase directory
    $RENAME_COMMAND = 'ren';		# DOS rename file command
    $XSUFX = '.exe';		# DOS suffix for executable module
    $PERLSLASH = $DOSSLASH;   # subdirectory delimitor in DOS file path
    $PLSUFX = '.pl';		# DOS suffix for perl script
    }
  elsif ( $UNkernel eq 'CYGWIN_NT' || $UNkernel eq 'MINGW32_NT' || $UNkernel eq 'MINGW64_NT' ) {
# unix-like environment under Windows
    $OStype = $UNkernel;     # set OStype for cygwin/mingw
# subdirectory delimitor for file paths is standard / under cygwin and
# under other unices
    $SLASH = $UNIXSLASH;      # subdirectory delimitor in file paths
    $PERLSLASH = $UNIXSLASH;    # subdirectory delimitor in file paths
    $COPY_COMMAND = 'cp -p';	# unix copy file command
    $DIFF_COMMAND = 'diff -b';	# unix difference command
    $ERASE_COMMAND = 'rm -f';	# unix delete file command
    $ERASE_RECURSIVE = 'rm -fr';	# unix command to recursively erase directories
    $RENAME_COMMAND = 'mv';	# unix rename file command
    $XSUFX = '';		# no special suffix for executable module
    $REDIRECT_ERR = '2>&1';	# unix shell gadget to redirect standard error
    $PLSUFX = '.pl';		# DOS suffix for perl script
    }
  else {
    print "unrecognized Windows Subsystem =$OStype=$UNkernel=";
    exit (99);
    }
  }
# if not Windows NT check other environment variables and uname output
elsif ($OStype eq 'OPENVMS') {
  $SLASH = $UNIXSLASH;                # subdirectory delimitor in file paths
  $COPY_COMMAND = 'copy';     # OpenVMS copy file command
  $ERASE_COMMAND = 'delete '; # OpenVMS delete file command
# for perl under normal unix systems:
  $XSUFX = '';		# use unix-style suffixes for binaries, ...
  $PLSUFX = '.pl';            # no special suffix for perl script under unix
  $PRLPFX = 'perl ';          # perl path defined in first line of script
  }
# MacOS section
elsif ($OStype eq 'MacOS') {
  $SLASH = $MACSLASH;             # subdirectory delimitor in file paths
  $COPY_COMMAND = 'Duplicate -c'; # copy file command
  $DIFF_COMMAND = 'Compare -b';   # file difference command
  $ERASE_COMMAND = 'Delete -y';   # delete file command
  $ERASE_RECURSIVE = 'Delete -y'; # command to recursively erase directories
  $SUFX = '';        		# MAC-style suffix for binaries
  $PERLSLASH = $MACSLASH;  	# subdirectory delimitor in MAC file paths
  $PLSUFX = '.pl';       	# suffix for perl script
  }
elsif ($OSname eq '' && $UNkernel ne '') {
# normal unix/linux section
  $OStype = $UNkernel;    	# set OStype for *ux
  $SLASH = $UNIXSLASH;          # subdirectory delimitor in file paths
  $COPY_COMMAND = '\cp -p';	# unix copy file command
  $DIFF_COMMAND = 'diff -b';	# unix difference command
  $ERASE_COMMAND = '\rm -f';	# unix delete file command
  $ERASE_RECURSIVE = '\rm -fr';	# unix command to recursively erase directories
  $RENAME_COMMAND = '\mv';	# unix rename file command
  $XSUFX = '';		# no special suffix for executable module
  $PERLSLASH = $UNIXSLASH;	# subdirectory delimitor for perl file paths
  $PLSUFX = '.pl';		# no special suffix for perl script
  }
else {
  print "unrecognized Operating System -$OSname=$OStype=$UNkernel-";
  exit (99);
  }
#
$CYGWIN = $UNkernel eq 'MINGW64_NT' ? '/cygwin' : '';
#
if ($TestDir ne '.') {
	if ( ! -e $TestDir )
	{
		print "Creating $TestDir";
		mkdir($TestDir);
	}
	print "cd $TestDir";
	chdir ("$TestDir") or die "could not access $TestDir";
	}
# open tests configuration file and print first line
$Config = "$top_testdir/$TestDir/tests.cnf";
$rc = open(CONF,"<$Config");
if ($rc eq '') {
	print "Error $rc opening file $Config";
	print 'Usage is: run-standard-tests machine_name serie first [last]';
	exit 20;
	}
$linect = 0;
while (<CONF>) {
	$linect ++;
	$X1 = substr($_,0,1);
	last if ($X1 ne '#');		# comments will be dropped
	}
print "Testing the Abinit code on the $MACH $OStype platform";
print $_;		# print title
print "Following tests will be run: $FIRST to $LAST";
#################################################
# set a date flag to make a new directory today #
#################################################
($sec,$min,$hour,$mday,$ymon,$yyear,$wday,$yday,$isdst)=localtime(time);
$ymon++;	# ymon was 0-11
$yyear +=1900;	# yyear was relative to 1900
$YYYYMMDD = sprintf("%4.4d",$yyear).sprintf("%2.2d",$ymon).sprintf("%2.2d",$mday);

$WORK_DIR = 'tmp-'.$MACH.'_'.$OStype.'_'.$YYYYMMDD;
if (! -e $WORK_DIR || ! -d $WORK_DIR) {
	mkdir ($WORK_DIR,0755);		# Mode 0755 ignored under DOS-Windows
	}
else {
	print "Do not create directory, $WORK_DIR already exists";
	}

print "cd $WORK_DIR";
chdir ("$WORK_DIR");

# import lib
do "$top_testdir/scripts/reportdiff.pl";

#Define fldiff.report file
$FLDREPORT = 'fldiff.report';	# report file for fldiff commands
# try to rename report file if it already exists
$rc = 0;
if (-f $FLDREPORT) {
	for ($i=1;$i <= 9;$i++) {
		$FLDREP2 = "old_fldiff$i".'.report';
		if (! -f $FLDREP2) {
			$rc = rename ($FLDREPORT,$FLDREP2);
			last;
			}
		}
	}
unlink ("$FLDREPORT") if ($rc == 0);	# get rid of file anyway

#Define chkinabi.report file
$CHKREPORT = 'tmp-chkinabi.report';   # report file for chkinabi commands
# try to rename report file if it already exists
$rc = 0;
if (-f $CHKREPORT) {
        for ($i=1;$i <= 9;$i++) {
                $CHKREP2 = "tmp-old_chkinabi$i".'.report';
                if (! -f $CHKREP2) {
                        $rc = rename ($CHKREPORT,$CHKREP2);
                        last;
                        }
                }
        }
unlink ("$CHKREPORT") if ($rc == 0);    # get rid of file anyway

#Define statrep file
if ( $TestDir eq "buildsys" )
{
	open(STATREPORT,">report");
}

# define file paths relative to $WORK_DIR in the unix fashion
# they will be translated by the perlpath/transpath function if necessary
$timeout_cmd = '';
#if ($timeout eq '0') {
#  $timeout_cmd = '';
#  }
#else {
#  $timeout_cmd = &perlpath("$abinit_nightlydir/timeout$XSUFX")." $timeout ";
#}
$CODE_SEQ = &perlpath("$abinit_bindir/abinit$XSUFX");	# codename for sequential version
$CODE_SEQ = "$timeout_cmd$CODE_SEQ";
$CODE_CPU = &perlpath("$abinit_bindir/abinit$XSUFX");  # sequential version for cpu tests
$CODE_CPU = "$timeout_cmd$CODE_CPU";
$CODE_MRGDDB = &perlpath("$abinit_bindir/mrgddb$XSUFX");	# codename for mrgddb
$CODE_MRGDDB = "$timeout_cmd$CODE_MRGDDB";
$CODE_MRGSCR = &perlpath("$abinit_bindir/mrgscr$XSUFX");	# codename for mrgscr
$CODE_MRGSCR = "$timeout_cmd$CODE_MRGSCR";
$CODE_MRGGKK = &perlpath("$abinit_bindir/mrggkk$XSUFX");	# codename for mrggkk
$CODE_MRGGKK = "$timeout_cmd$CODE_MRGGKK";
$CODE_ANADDB = &perlpath("$abinit_bindir/anaddb$XSUFX");	# codename for anaddb
$CODE_ANADDB = "$timeout_cmd$CODE_ANADDB";
$CODE_CUT3D = &perlpath("$abinit_bindir/cut3d$XSUFX");	# codename for cut3d
$CODE_CUT3D = "$timeout_cmd$CODE_CUT3D";
$CODE_AIM = &perlpath("$abinit_bindir/aim$XSUFX");	# codename for atom-in-molecule (aim) code
$CODE_AIM = "$timeout_cmd$CODE_AIM";
$CODE_CONDUCTI = &perlpath("$abinit_bindir/conducti$XSUFX");	# codename for conducti
$CODE_CONDUCTI = "$timeout_cmd$CODE_CONDUCTI";
$CODE_OPTIC = &perlpath("$abinit_bindir/optic$XSUFX");	# codename for optic
$CODE_OPTIC = "$timeout_cmd$CODE_OPTIC";
$CODE_MACROAVE = &perlpath("$abinit_bindir/macroave$XSUFX");	# codename for macroave
$CODE_MACROAVE = "$timeout_cmd$CODE_MACROAVE";
$CODE_UJDET = &perlpath("$abinit_bindir/ujdet$XSUFX");    # codename for ujdet
$CODE_UJDET = "$timeout_cmd$CODE_UJDET";
$CODE_LWF = &perlpath("$abinit_bindir/lwf$XSUFX");	# codename for lwf
$CODE_LWF = "$timeout_cmd$CODE_LWF";
$CODE_BAND2EPS = &perlpath("$abinit_bindir/band2eps$XSUFX");	# codename for band2eps
$CODE_BAND2EPS = "$timeout_cmd$CODE_BAND2EPS";
$CODE_ATOMPAW = &perlpath($ENV{"run_atompaw"});    # codename for atompaw
$CODE_ATOMPAW = "$timeout_cmd$CODE_ATOMPAW";
$CODE_VDWKG = &perlpath("$abinit_bindir/vdw_kernelgen$XSUFX");	# codename for vdw_kernelgen
$CODE_VDWKG = "$timeout_cmd$CODE_VDWKG";
$CODE_DOCCHK = &perlpath("$top_testdir/scripts/docchk.py$XSUFX");    # codename for docchk
$CODE_DOCCHK = "$timeout_cmd$CODE_DOCCHK";
$CODE_GUI = &perlpath("$java_runner -jar $top_srcdir/../gui/precompiled/AbinitGUI.jar$XSUFX");    # codename for gui
$CODE_GUI = "$timeout_cmd$CODE_GUI";
$CODE_CHECK_UNALLOWED = &perlpath("$top_testdir/scripts/check_forbidden.py$XSUFX");    # codename for check_forbidden
$CODE_CHECK_UNALLOWED = "$timeout_cmd$CODE_CHECK_UNALLOWED";
$CODE_CHECK_ASCII = &perlpath("$top_testdir/scripts/check_ascii.py$XSUFX");  # codename for check_ascii
$CODE_CHECK_ASCII = "$timeout_cmd$CODE_CHECK_ASCII";
$CODE_CHECK_INLINED_MACROS = &perlpath("$top_testdir/scripts/check_inlined_macros.py$XSUFX");  # codename for check_inlined_macros
$CODE_CHECK_INLINED_MACROS = "$timeout_cmd$CODE_CHECK_INLINED_MACROS";
$CODE_CHECK_FORBIDDEN_IN_DOC_YML = &perlpath("$top_testdir/scripts/check_forbidden_in_doc_yml.py$XSUFX");  # codename for check_forbidden_in_doc_yml
$CODE_CHECK_FORBIDDEN_IN_DOC_YML = "$timeout_cmd$CODE_CHECK_FORBIDDEN_IN_DOC_YML";
$CODE_WARNCHK = &perlpath("$top_testdir/scripts/warningschk.py$XSUFX");    # codename for warnchk
$CODE_WARNCHK = "$timeout_cmd$CODE_WARNCHK";
$CODE_CHKINPVARS = &perlpath("$top_testdir/scripts/check-input-vars.py$XSUFX");    # codename for checkinputvars
$CODE_CHKINPVARS = "$timeout_cmd$CODE_CHKINPVARS";
$CODE_FFTPROF = &perlpath("$abinit_bindir/fftprof$XSUFX");    # codename for fftprof
$CODE_FFTPROF = "$timeout_cmd$CODE_FFTPROF";
$CODE_PAR = '/usr/local/mpi-pgi4/bin/mpirun -np 2 -machinefile sleepy.pcpm.ucl.ac.be:2 ../../../src/98_main/abinit' ;      # run parallel version on sleepy
$CODE_PAR = "$timeout_cmd$CODE_PAR";

$SORTED_CMD = &perlpath("$PERL $top_testdir/scripts/Sort.sh");	# sorting lines for warningschk.py : t03 & t04

$INPUTDIR = &transpath("$CYGWIN$top_testdir/$TestDir/Input");	# input directory
$REF = &transpath("$top_testdir/$TestDir/Refs");		# reference directory
$PSPS = &transpath("$CYGWIN$top_testdir/Psps_for_tests");	# pseudopotential directory
$CHKINABI = &perlpath("$PERL $top_testdir/scripts/chkinabi.pl");  # relative paths
$FLDIFF = &perlpath("$PERL $top_testdir/scripts/fldiff.pl");	# to scripts directory
#
# Read tests configuration file until first test is found. The file format is:
#			TestID [command|program] fn1=path1 fn2=path2 ... [opt=Opt1] ...
# TestID is a 1 to 3 characters string
# command may be "copy" or "erase"
# program refers to an ABINIT program, valid names are:
# 	abinit | seq	exercise sequential version
#		aim
#		atompaw
#		ana | anaddb
#		band2eps
#		chi
#		conducti
#		fftprof
#		optic
#		cut3d
#		lwf
#		macroave
#		mrg | mrgddb
#		mrggkk
#		mrgscr
#               ncdump
#               ujdet
#
#   For all the tests, use unix-style file paths, they will be converted
# by function transpath if necessary. File paths can contain the perl variable
# $REF or $PSPS.
#   See the routine setfnopt for valid options
#

$CurTest = '';
while (<CONF>) {
	$linect ++;
	$X1 = substr($_,0,1);
	next if ($X1 eq '#');		# comments will be dropped
	($testID,$prgcmd,$P1,$P2,$P3,$P4,$P5,$P6,$P7,$P8,$P9,$P10,$P11,$P12) = split(' ',$_);
	next if ($testID eq '' || length($testID) > 10 || $prgcmd eq '');	# handle as comments
	if ($CurTest eq '') {
		if ($testID eq $FIRST || $FIRST eq 'first') {
			$CurTest = $testID;	# first test has been found
			$option = 0;		# option for abinit tests
			$fldopt = '';			# extra options for fldiff script, default is none
			$fldflag = $fldRun;		# diff/fldiff default option
			}
		else {
			next;		# skip until first test
			}
		}
	if ($testID ne $CurTest){		# first line of another test ?
		last if ($CurTest eq $LAST || ($testID eq 'end' && $LAST eq 'end'));	# all done
		$CurTest = $testID;		# prepare for next test
		$option = 0;		# option for abinit tests
		$fldopt = '';			# extra options for fldiff script, default is none
		$fldflag = $fldRun;		# diff/fldiff default option
		}
# ******************************************************************
# see the comments of each section for explanations of the operations                                                   #
# ******************************************************************
	print "Test $CurTest pgm/cmd $prgcmd $P1 $P2 $P3 ..." if ($debug >= 2);
	if ($prgcmd eq 'seq' || $prgcmd eq 'abinit') {
		&doseqtest($CurTest,$P1,$P2,$P3,$P4,$P5,$P6,$P7,$P8,$P9,$P10,$P11,$P12);
		}
	elsif ($prgcmd eq 'ana' || $prgcmd eq 'anaddb') {
		&doanatest($CurTest,$P1,$P2,$P3,$P4);
		}
	elsif ($prgcmd eq 'mrg' || $prgcmd eq 'mrgddb') {
		$ix = index($_,$P1);
		$ParmString = substr($_,$ix);
		&domrgtest($CurTest,$ParmString);
		}
        elsif ($prgcmd eq 'mrgscr') {
                &domrgscrtest($CurTest);
		}
        elsif ($prgcmd eq 'mrggkk') {
                &domrggkktest($CurTest);
                }
	elsif ($prgcmd eq 'cut3d') {
		&docut3dtest($CurTest,$P1,$P2,$P3,$P4,$P5);
		}
        elsif ($prgcmd eq 'ujdet') {
                &doujdettest($CurTest,$P1,$P2,$P3,$P4,$P5);
                }
	elsif ($prgcmd eq 'aim') {
		&doaimtest($CurTest,$P1,$P2,$P3,$P4,$P5,$P6);
		}
	elsif ($prgcmd eq 'conducti') {
		&doconductitest($CurTest,$P1,$P2,$P3);
		}
        elsif ($prgcmd eq 'fftprof') {
                &dofftproftest($CurTest);
                }
        elsif ($prgcmd eq 'optic') {
                &dooptictest($CurTest,$P1,$P2,$P3);
                }
	elsif ($prgcmd eq 'lwf') {
		&dolwftest($CurTest,$P1,$P2,$P3,$P4,$P5);
		}
        elsif ($prgcmd eq 'atompaw') {
                &doatompawtest($CurTest,$P1);
                }
        elsif ($prgcmd eq 'gui') {
                &doguitest($CurTest,$P1);
                }
        elsif ($prgcmd eq 'vdw_kgen') {
                &dovdwkgtest($CurTest,$P1);
                }
        elsif ($prgcmd eq 'docchk') {
                &dodocchk($CurTest);
                }
        elsif ($prgcmd eq 'check_forbidden') {
                &docheck_forbidden($CurTest);
                }
        elsif ($prgcmd eq 'check_ascii') {
                &docheck_ascii($CurTest);
                }
        elsif ($prgcmd eq 'check_inlined_macros') {
                &docheck_inlined_macros($CurTest);
                }
        elsif ($prgcmd eq 'check_forbidden_in_doc_yml') {
                &docheck_forbidden_in_doc_yml($CurTest);
                }
        elsif ($prgcmd eq 'chkinpvars') {
                &dochkinpvars($CurTest);
                }
        elsif ($prgcmd eq 'warnchk') {
                &dochkwarnings($CurTest,$P1);
                }
	elsif ($prgcmd eq 'band2eps') {
		&doband2epstest($CurTest,$P1,$P2,$P3,$P4,$P5);
		}
	elsif ($prgcmd eq 'macroave') {
		&domacroavetest($CurTest,$P1,$P2,$P3,$P4,$P5);
		}
	elsif ($prgcmd eq 'cpuserie') {
		&doserieX($CurTest,$P1,$P2,$P3,$P4,$P5);
		}
	elsif ($prgcmd eq 'copy') {
	  $P1 =~ s/\$PSPS/$PSPS/g;
	  $P1 =~ s/\$INPUTDIR/$INPUTDIR/g;
	  $P1 =~ s/\$REF/$REF/g;
		&docopy($CurTest,$P1,$P2);
		}
	elsif ($prgcmd eq 'erase') {
		&doerase($CurTest,$P1);
		}
        elsif ($prgcmd eq 'report') {
                &doreport($CurTest);
                }
        elsif ($prgcmd eq 'ncdump') {
                &doncd($CurTest,$P1,$P2);
                }
        elsif ($prgcmd eq 'statchk') {
                &dochkstatus($CurTest,$P1);
                }

# end of while loop
	else {
		print "Unknown program $prgcmd found at line $linect with test $CurTest"
		}
	}
if ($CurTest eq '') {		# end of file ?
	print "Error: test $FIRST is not defined.";
	exit 8;
	}
End_of_tests:
#Close statrep file
if ( $TestDir eq "buildsys" )
{
	close(STATREPORT);
}
exit 0;

# ****************************************
sub transpath {
	local($path) = @_;
#
# purpose: translate unix-like path of file to DOS-like according to $SLASH
# argument:
#	$path = path to be translated
# output: this subroutine acts as a function and returns the path
# according to host conventions

	$path =~ tr/\057/\\/ if ($SLASH eq '\\');
	if ($OStype eq 'MacOS') {
		$path =~ tr/\057/:/;
		$path =~ s/\.\.//g;
		$path =~ s/:/::/ ;
		}
	return $path;
	}
# ****************************************
sub perlpath {
	local($path) = @_;
#
# purpose: translate unix-like path of file to DOS-like according to $PERLSLASH.
#   This is necessary when calling a DOS command like perl under PGI Workstation.
# argument:
#	$path = path to be translated
# output: this subroutine acts as a function and returns the path
# according to host conventions

	$path =~ tr/\057/\\/ if ($PERLSLASH eq '\\');
	if ($OStype eq 'MacOS') {
		$path =~ tr/\057/:/;
		$path =~ s/\.\.//g;
		$path =~ s/:/::/ ;
		}

	return $path;
	}
# ****************************************
sub setfnopt {
	local($parm) = @_;
#
#  purpose: decode the argument and set the corresponding variable
# argument: format is "fnI=path" for the file names, "pspI=path" for
# pseudo-potential files, "GEO=path" for the GEO file, "DS3AT=path"
# for the DS3_DOS_AT000x file fftalg=abc for the cpu tests FFT algorithm
# or files=fn1,fn2[,...,fnLast] for a mrgddb or mrgscr files-to-be-merged list
# or "opt=OPTION" for options
#  output: either the variable $fnI, $psp[I], $GEOfn is set to "path"
# or $FFTalg, or $fldopt, or $fldflag, or $option is set
#  options:
# a) standard options of the fldiff script (see the script source)
# will be passed to it when called:
#		-easy		-medium	 -ridiculous
# b) options to tailor the analysis of the results after test completion
#		-nofld	-diflog	-DOS	-DS2	+xml
# c) other options to modify the flow of operations during the test
#		-keepfn4 -istwfk1 -istwfk11
	$_ = $parm;
# try fnI=
	$hit1 = m/fn([1-9])=(.*)/;
	if ($hit1) {
# build perl instruction to assign new value to $fnI
		$isn = '$fn'.$1.'="'.$2.'"';
		print $isn if ($debug >= 3);
# assign new file name to $fnI
# substitution of perl variables like $REF ou $PSPS will be done
		eval($isn);
		return;
		}
# try pspI=
	$hit2 = m/psp([1-9])=(.*)/;
	if ($hit2 && ($prgcmd eq 'seq' || $prgcmd eq 'abinit' || $prgcmd eq 'aim')) {
# assign new pseudo-potential file name to $psp[i]
		$psp[$1] = "$PSPS/$2";
		return;
		}
# try pspatompawI=
        $hit2 = m/pspatompaw([1-9])=(.*)/;
        if ($hit2 && ($prgcmd eq 'seq' || $prgcmd eq 'abinit' || $prgcmd eq 'aim')) {
# assign new pseudo-potential file name to $psp[i]
                $psp[$1] = "$2";
                return;
                }
# try GEO=
	$hit3 = m/GEO=(.*)/;
	if ($hit3 && ($prgcmd eq 'seq' || $prgcmd eq 'abinit')) {
# assign new file name to GEO file
		$GEOfn = $1;
		return;
		}
# try DS3AT=
	$hit3 = m/DS3AT=(.*)/;
	if ($hit3 && ($prgcmd eq 'seq' || $prgcmd eq 'abinit')) {
# assign new file name to DS3_DOS_AT000x file
		$DS3ATfn = $1;
		return;
		}
# try files=
	$hit3 = m/files=(.*)/;
	if ($hit3 && ($prgcmd eq 'mrgddb' || $prgcmd eq 'mrg' || $prgcmd eq 'mrgscr')) {
		@mrgfiles = split(',',$1);
		}
# try fftalg=
	$hit4 = m/fftalg=(.*)/;
	if ($hit4 && $prgcmd eq 'cpuserie') {
# set new value for FFT algorithm
		$FFTalg = $1;
		return;
		}
# try option, 1st for cut3d
	$hit5 = m/opt=(.*)/;
	if ($hit5 && $1 eq '-logout' && $prgcmd eq 'cut3d') {
		$fldflag = $fldlog;	# run diff/fldiff on the log file
		return;
		}
        elsif ($1 eq '-nofld' && $prgcmd eq 'cut3d') {
                        $fldflag = $fldNorun;   # don't run the fldiff script
                }

# try option for cpu tests
	if ($hit5 && $prgcmd eq 'cpuserie') {
		if ($1 eq '-istwfk1') {
			$option = $optistwfk1;
			}
		elsif ($1 eq '-istwfk11') {
			$option = $optistwfk11;
			}
		else {
			print "Error, invalid parameter $parm at line $linect with test $CurTest";
			}
		return;
		}
# try fldiff option for abinit, aim, cut3d, mrgddb, mrgscr or lwf
	return if ($prgcmd ne 'seq' && $prgcmd ne 'abinit' && $prgcmd ne 'cut3d' && $prgcmd ne 'aim' &&
		$prgcmd ne 'lwf' && $prgcmd ne 'macroave' && $prgcmd ne 'anaddb' && $prgcmd ne 'conducti' && $prgcmd ne 'fftprof' &&
		$prgcmd ne 'mrggkk' && $prgcmd ne 'optic' && $prgcmd ne 'mrgddb' && $prgcmd ne 'mrg' && $prgcmd ne 'mrgscr');
	if ($hit5 && ($1 eq '-easy' || $1 eq '-medium' || $1 eq '-ridiculous') ) {
		$fldopt = $1;
		return;
		}
	return if ($prgcmd ne 'seq' && $prgcmd ne 'abinit');
# try abinit option
	if ($hit5) {
		if ($1 eq '-nofld') {
			$fldflag = $fldNorun;	# don't run the fldiff script
			}
		elsif ($1 eq '-DOS') {
			$fldflag = $fldDOS;	# run diff/fldiff on the o_DOS file
			}
		elsif ($1 eq '-DS2') {
			$fldflag = $fldDS2;	# run diff/fldiff on the o_DS2_DOS file
			}
		elsif ($1 eq '+xml') {
			$fldflag = $fldxml;	# run diff/fldiff also on the xml file
			}
		elsif ($1 eq '-keepfn4') {
			$option = $optKeepfn4;	# keep the output wfs files (fn4*) after completion
			}
		elsif ($1 eq '-logxml') {
			$fldflag = $fldlgx;	# run diff/fldiff on the o_LOG.xml file
			}
		elsif ($1 eq '-MDF') {
			$fldflag = $fldMDF;	# run diff/fldiff on the MDF files
			}
		else {
			print "Error, invalid option $1 at line $linect with test $CurTest";
			}
		return;
		}
	print "Error, invalid parameter $parm at line $linect with test $CurTest";
	}
# ****************************************
sub copyunzip {
	local($path) = @_;
#
# purpose: copy file from parent to current directory and g-unzip it if .gz-suffixed
# argument:
#	$path = path of file
# output: this subroutine acts as a function and returns the unzipped file name
# Remove the $CYGWIN prefix for system commands like copy diff...
        if ( $path =~ m/\/cygwin(.*)/ ) { $path = $1; };
        $path = $1 if ($& ne '');
#
	$path2 = &transpath("$path");	# translate path name when required
	$path3 = $path;
# if the file doesn't exist, try the g-zipped file
	if (! -f $path2) {
		$path2 = $path2.'.gz';
		$path3 = $path.'.gz';
		}
	$rc = system("$COPY_COMMAND $path2 .");	# copy file in working directory
	if ($rc != 0) {
		print "Error $rc copying file $path";
		return '';
		}
	$rix = rindex($path3,$UNIXSLASH);	# find last subdirectory delimitor
	$path2 = $rix <= 0 ? $path3 : substr($path3,$rix+1);	# extract file name
	$len3 = length($path2) - 3;
	$sufx = substr($path2,$len3,3);	# pick up last 3 characters
	if ($sufx eq '.gz') {		# g-zipped file ?
# Since .gz format is uncommon on some non-unix systems (e.g. Win9x), the following
# command may fail:
		$rc = system("gzip -d $path2");		# g-unzip file
		if ($rc != 0) {
			print "Error $rc unzipping file $path2";
			return '';
			}
		else {
			print "File $path copied to $path2 and unzipped";
			return substr($path2,0,$len3);	# return unzipped file name
			}
		}
	else {
		print "File $path copied to $path2";
		}
	return $path2;		# return file name
	}
# ****************************************
sub dofldiff {
	local($label,$result,$reference,$fldopt) = @_;
#
# purpose: call the fldiff perl script for a floating point comparison of the result
#	file with reference file and write output to $FLDREPORT (fldiff.report).
# arguments:
#	$label = label to be written to $FLDREPORT, usually "Case_$TID"	# Case_01
#	$result = test output file path
#	$reference = reference file path
# $fldopt = specific fldiff options for the Case
	return if ($OStype eq 'MacOS');
# append label with test number to report file
	open (FLDREP,">>$FLDREPORT") || die "Unable to open FLDREP for $label";
	print FLDREP "\n","$label";
	close (FLDREP);
# Use the floating diff script to get a compact report on the Run
	print "Doing floating point compare of $result with $reference and option $fldopt";
	$result = &transpath($result);	# translate path name according to $PERLSLASH
	$reference = &transpath($reference);
	print "$FLDIFF $deffldopt $fldopt $result $reference $label" if ($debug >= 2); # GAF
	system ("$FLDIFF $deffldopt $fldopt $result $reference $label >> $FLDREPORT 2>> report.in");
	return;
	}
# ****************************************
sub doseqtest {
	local($TID,$p1,$p2,$p3,$p4,$p5,$p6,$p7,$p8,$p9,$p10,$p11,$p12) = @_;
#
# purpose: exercise the sequential version
# arguments:
#		$TID = test identification (1 to 3 characters string)
#		$p1, ... = 1 to 12 parameters [re]defining some options or the file names
# (fn1 to fn9) that will be put into the "files" file named "abfiles.run"
# this file will be used as input to abinit
#	fn2 is the output file name that will be compared with reference
#	fn5 will be used as suffix to diff and prefix to log file names
# pseudo-potential file names can be specified for fn6-fn9, as psp1-psp4
# the parameter format for a file name is: fnI=path
# the parameter format for a pseudo-potential file name is: pspI=path
# the parameter format for the GEO file name is: GEO=path
# the parameter format for the DS3_DOS_AT000x file name is: DS3AT=path
# the parameter format for an option is: opt=option ; see setfnopt
# set default files names derived from $TID:
	$fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";	# e.g. ../t01.in
	$fn2 = "t$TID.out";		# t01.out
	$fn3 = "t$TID".'i';		# t01i
	$fn4 = "t$TID".'o';		# t01o
	$fn5 = "t$TID";				# t01
	$fn6 = '';		# undefined
	$fn7 = '';		# undefined
	$fn8 = '';		# undefined
	$fn9 = '';		# undefined
	$psp[1] = '';		# undefined
	$psp[2] = '';		# undefined
	$psp[3] = '';		# undefined
	$psp[4] = '';		# undefined
	$GEOfn = '';		# undefined
	$DS3ATfn = '';		# undefined
#DEBUG
#        print "$fn5";
# set new file names or options according to the parameters
	for ($i=1;$i <= 9;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#DEBUG
#        print "$fn5";
	$fn6 = $psp[1] if ($psp[1] ne '');
	$fn7 = $psp[2] if ($psp[2] ne '');
	$fn8 = $psp[3] if ($psp[3] ne '');
	$fn9 = $psp[4] if ($psp[4] ne '');
  print "fldflag $fldflag fldopt $fldopt" if ($debug >= 3);
#	print "Case_$TID:";		# Case_01
#DEBUG
#        print "$fn5";
	print $fn1,$fn2,$fn3,$fn4,$fn5,$fn6,$fn7,$fn8,$fn9 if ($debug >= 2);
# Test input file in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                system ("$CHKINABI $fn1 >> $CHKREPORT");
                return ;
                }
#DEBUG
#        print "$fn5";
# Remove existing files related to this test case:
# In case of --without-stdout option, we erase the "log" file also.
	$RUNfile = 'ab.files';		# files file for abinit
	$logfn = $fn5.'.log';		# t01.log
	$errfn = $fn5.'.err';		# t01.err
	$difffn = 'diff.'.$fn5;		# diff.t01
	$difibmfn = 'difibm.'.$TID;	# difibm.01
	$outfn = &transpath($fn2);	# translate path name when required
	unlink($RUNfile,$difffn,$difibmfn,$logfn,$errfn,$outfn,"log");
# For v1 Case 63 the output wfs file should not be erased
	if ($option != $optKeepfn4) {
		$fnglob = $fn4."_$MULTI";		# t01o_*
		system("$ERASE_COMMAND $fnglob");	# unlink <t01o_*>
		}
	$fnglob = $fn5."_$MULTI";		# t01_*
	system("$ERASE_COMMAND $fnglob");	# unlink < t01_* >

# Create "files" file with file names $fn1, ... one per line
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES &transpath($fn1);	# translate path name when required
	print FILES $outfn;
	print FILES $fn3;
	print FILES $fn4;
	print FILES $fn5;
	print FILES &transpath($fn6);	# translate path name when required
	print FILES &transpath($fn7) if ($fn7 ne '');	# add 7th filename if any
	print FILES &transpath($fn8) if ($fn8 ne '');	# add 8th filename if any
	print FILES &transpath($fn9) if ($fn9 ne '');	# add 9th filename if any
	close(FILES);
# Run the test case
	return if ($debug >= 3);
	print "\n[$TestDir][$fn5] File $RUNfile created, starting abinit";
	$Time_Before=time();
	if ($OStype eq 'MacOS') {
		system ("$CODE_SEQ < $RUNfile \xb7 $logfn ");
		}
	elsif ($MACH eq 'paral') {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_PAR < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_SEQ < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
	$Time_After=time();
	$Runtime=$Time_After-$Time_Before;
# Compare with reference files
	$REFoutfn = &transpath("$REF/$fn2");
	print "[$TestDir][$fn5] Finished abinit (Runtime: $Runtime seconds)";
        print "[$TestDir][$fn5] Comparing $outfn with reference file";
	print "[$TestDir][$fn5] Reference file: $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
# handle o_DOS file for v1 case #07
	if ($fldflag == $fldDOS) {
		$oDOSfile = "t$TID".'o_DOS';	# t07o_DOS
		$REFoDOS = &transpath("$REF/$oDOSfile");
		print "[$TestDir][$fn5] Comparing $oDOSfile with $REFoDOS";
		$difoDOS = $difffn.'_DOS';		# diff.t07_DOS
		system ("$DIFF_COMMAND $oDOSfile $REFoDOS > $difoDOS");
		&dofldiff("Case_$TID : DOS files",$oDOSfile,$REFoDOS,$fldopt);	# compare floating point numbers
		}
# handle o_DS2_DOS file for v3 case #46, v4 case #35
	elsif ($fldflag == $fldDS2) {
		$oDOSfile = "t$TID".'o_DS2_DOS';	# t46o_DS2_DOS
		$REFoDOS = &transpath("$REF/$oDOSfile");
		print "[$TestDir][$fn5] Comparing $oDOSfile with $REFoDOS";
		$difoDOS = $difffn.'_DS2_DOS';              # diff.t46_DS2_DOS
		system ("$DIFF_COMMAND $oDOSfile $REFoDOS > $difoDOS");
		&dofldiff("Case_$TID : DOS file",$oDOSfile,$REFoDOS,$fldopt);
		}
	elsif ($fldflag == $fldlgx) {
# Use the floating diff script on the _LOG.xml file, used in v5
		$oXMLfile = "t$TID".'o_LOG.xml';
		$REFoXML = &transpath("$REF/$oXMLfile");
		print "Comparing $oXMLfile with $REFoXML";
		$difoXML = $difffn.'_LOG';
		system ("$DIFF_COMMAND $oXMLfile $REFoXML > $difoXML");
		&dofldiff("Case_$TID : LOG.xml file",$oXMLfile,$REFoXML,$fldopt);
		}
	elsif ($fldflag == $fldMDF) {
# Use the floating diff script on all the MDF files produced by the test.
                #
                # Compare the main output file.
	        &dofldiff("Case_$TID :",$outfn,$REFoutfn,$fldopt);
                #
                # Compare the MDF files.
                $pattern = "t$TID*_MDF";
                my @mdf_files = glob($pattern);

                foreach my $out_mdf (@mdf_files) {
		  $ref_mdf = &transpath("$REF/$out_mdf");
                  print "Comparing MDF file $out_mdf with reference $ref_mdf\n";
		  $difomdf = "diff_".$out_mdf;
		  system ("$DIFF_COMMAND $out_mdf $ref_mdf > $difomdf");
		  &dofldiff("Case_$TID : $out_mdf vs $ref_mdf",$out_mdf,$ref_mdf,$fldopt);
                  }
                }
# v1 Case 98 and 99: no use of fldiff for these tests
	elsif ($fldflag != $fldNorun) {
# Use the floating diff script to get a compact report on the Run
		$REFoutfn = "$REF/$fn2";
		&dofldiff("Case_$TID",$fn2,$REFoutfn,$fldopt);	# compare floating point numbers
		}
# handle DS3_DOS_AT000x file for tests v4 #35,38
	if ($DS3ATfn ne '') {
		$REFDS3AT = &transpath("$REF/$DS3ATfn");	# t35o_DS3_DOS_AT0001
		$difoDOS = $difffn.substr($DS3ATfn,4);		# diff.t35_DS3_DOS_AT0001
		print "[$TestDir][$fn5] Comparing $DS3ATfn with $REFDS3AT";
		system ("$DIFF_COMMAND $DS3ATfn $REFDS3AT > $difoDOS");
		&dofldiff("Case_$TID : DOS file",$DS3ATfn,$REFDS3AT,$fldopt);	# compare floating point numbers
		}
# handle GEO file for tests fast #28, 29: compare with reference
	if ($GEOfn ne '') {
# compare with reference file
		$diffn = 'diff.'.$TID.'_GEO';		# diff.28_GEO
		$REFGEOfn = &transpath("$REF/$GEOfn");
		print "[$TestDir][$fn5] Comparing $GEOfn with $REFGEOfn";
		system ("$DIFF_COMMAND $GEOfn $REFGEOfn > $diffn");
		&dofldiff("Case_$TID : GEO file",$GEOfn,$REFGEOfn,$fldopt);	# compare floating point numbers
		}
# v1 Case 40, 41 and 42: also take into account the xml files
	if ($fldflag == $fldxml) {
		$cmlfn = $fn4.'.cml';	# t40o.cml
		$REFcmlfn = &transpath("$REF/$cmlfn");
		print "[$TestDir][$fn5] Comparing $cmlfn with $REFcmlfn";
		$difffn = 'diff.'.$cmlfn;         # diff.t40o.cml
		unlink($difffn);
		system ("$DIFF_COMMAND $cmlfn $REFcmlfn > $difffn");
		&dofldiff("Case_$TID : .cml file",$cmlfn,$REFcmlfn,$fldopt);	# compare floating point numbers
		}
	return;			# go read next line from configuration file
	}
# ****************************************
sub doanatest {
	local($TID,$p1,$p2,$p3,$p4) = @_;
	local($fntodelete);
#
# purpose: exercise anaddb
# arguments:
#		$TID = test identification (1 to 3 characters string)
#	that will be used as suffix to diff and prefix to log file names
#		$p1, ...,$p4 = 1 to 4 parameters [re]defining the file names
#   fn1 to fn6 that will be put into the "files" file named "anaddb.run";
# this file will be used as input to anaddb;
# the parameters format is: fnI=path
#   fn2 is also the output file name that will be compared with reference.
#   fn3 can reside in the Test directory as a normal or g-zipped file,
# it will be copied, and unzipped if necessary, to the working directory before
# the test and the copy will be deleted after
#   set default file name(s) derived from $TID:
	$fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";	# e.g. ../t13.in
	$fn2 = "t$TID.out";	# t13.out
	$fn3 = "t$TID.ddb.in";	# t13.ddb.in
	$fn4 = "t$TID.md";	# t13.md
	$fn5 = "t$TID.gkk";	# t13.gkk Dummy file for now
	$fn6 = "t$TID";		# t13  base name for elphon outputs
	$fn7 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.ddk";	# ../t13.ddk default name for file containing list of ddk filenames
# set new file names according to the parameters
# Added mverstra 29 07 2004
	for ($i=1;$i <= 5;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#
#	print "Case_$TID:";
# Test input file in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                system ("$CHKINABI -m anaddb $fn1 >> $CHKREPORT");
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'anaddb.run';		# files file for anaddb
	$logfn = 't'.$TID.'.log';		# t13.log
	$errfn = 't'.$TID.'.err';		# t13.log
	$difffn = 'diff.t'.$TID;		# diff.t13
	$outfn = &transpath($fn2);	# translate path name when required
	unlink($RUNfile,$difffn,$logfn,$errfn,$outfn);
# if file fn3 does not exist, it is probably to be copied/unzipped from
# the parent directory
	if (! -f $fn3) {
		$fntodelete = &copyunzip("$top_testdir/$TestDir/Input/$fn3");	# copy/unzip ../t13.ddb.in[.gz]
		}
	else {
		$fntodelete = '';	# undefined
		}
# Create "files" file with file names $fn1, ... one per line
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES &transpath($fn1);	# translate path name when required
	print FILES $outfn;
	print FILES $fn3;
	print FILES $fn4;
	print FILES $fn5;
	print FILES $fn6;
	print FILES &transpath($fn7);
	close(FILES);

	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] File $RUNfile created, starting anaddb";
	if ($OStype eq 'MacOS') {
		system ("$CODE_ANADDB < $RUNfile \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_ANADDB < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
# Compare with reference output file
	$REFoutfn = &transpath("$REF/$fn2");
	print "[$TestDir][t$TID] Finished anaddb, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);	# compare floating point numbers
# Compare with reference log file
#       if("$TestDir/$TID" eq 'v5/28') {
#       $REFlogfn = &transpath("$REF/$logfn");
#       print "Finished anaddb, comparing $logfn with $REFlogfn";
#       system ("$DIFF_COMMAND $logfn $REFlogfn > 'diff.t28.log'");
#       &dofldiff("Case_$TID.log","$logfn","$REF/$logfn",$fldopt);      # compare floating point numbers
#             }
# check for work file to be deleted
	system("$ERASE_COMMAND $fntodelete") if ($fntodelete ne '');	# delete file
	return;			# go read next line from configuration file
	}

# ****************************************
sub doatompawtest {
        local($TID) = @_;
#
# purpose: exercise atompaw
# arguments:
#               $TID = test identification that will be used as prefix to input file names
# the parameter format for a file name is: fnI=path
# $fn1 = input to atompaw
# set default files names derived from $TID:
        $fn1 = "$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../telphon_3.in
        $infn = &transpath($fn1);
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $logfn = 't'.$TID.'.log';               # t13.log
        $errfn = 't'.$TID.'.err';               # t13.log
        unlink($logfn,$errfn);

# Run the test case
        print "\n[$TestDir][t$TID] Starting atompaw";
        if ($OStype eq 'MacOS') {
                system ("$CODE_ATOMPAW < $infn \xb7 $logfn ");
                }
        else {
                $REDIRECT_ERR = "2> $errfn";
                system ("$CODE_ATOMPAW < $infn > $logfn $REDIRECT_ERR");
                $REDIRECT_ERR = '2>&1';
                }
        print "[$TestDir][t$TID] Finished atompaw";
        return;                 # go read next line from configuration file
        }

# ****************************************
sub dovdwkgtest {
        local($TID) = @_;
#
# purpose: exercise atompaw
# arguments:
#               $TID = test identification that will be used as prefix to input file names
# the parameter format for a file name is: fnI=path
# $fn1 = input to atompaw
# set default files names derived from $TID:
        $fn1 = "$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../telphon_3.in
        $infn = &transpath($fn1);
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = 't'.$TID.'.out';               # t13.out
        $logfn = 't'.$TID.'.log';               # t13.log
        $errfn = 't'.$TID.'.err';               # t13.err
	$difffn = "diff.t$TID";			# diff.t13
        unlink($logfn,$errfn);

# Run the test case
        print "\n[$TestDir][t$TID] Starting vdw_kernelgen";
        if ($OStype eq 'MacOS') {
                system ("$CODE_VDWKG < $infn \xb7 $logfn ");
                }
        else {  
                $REDIRECT_ERR = "2> $errfn";
                system ("$CODE_VDWKG < $infn > $logfn $REDIRECT_ERR");
                $REDIRECT_ERR = '2>&1';
                }
        system("cp $logfn $outfn");
# Compare with reference files
	$REFoutfn = &transpath("$REF/$outfn");
	print "[$TestDir][t$TID] Finished vdw_kernelgen, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$outfn","$REFoutfn",$fldopt);	# compare floating point numbers
        return;                 # go read next line from configuration file
        }

# ****************************************
sub doguitest {
        local($TID) = @_;
#
# purpose: exercise gui
# arguments:
#               $TID = test identification that will be used as prefix to input file names
# the parameter format for a file name is: fnI=path
# $fn1 = input to atompaw
# set default files names derived from $TID:
        $fn1 = "$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../t01.in
        $infn = &transpath($fn1);
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $logfn = 't'.$TID.'.log';               # t01.log
        $errfn = 't'.$TID.'.err';               # t01.log
        unlink($logfn,$errfn);

# Run the test case
        print "\n[$TestDir][t$TID] Starting GUI";
        if ($OStype eq 'MacOS') {
                system ("$CODE_GUI < $infn \xb7 $logfn ");
                }
        else {
                $REDIRECT_ERR = "2> $errfn";
                system ("$CODE_GUI < $infn > $logfn $REDIRECT_ERR");
                $REDIRECT_ERR = '2>&1';
                }
        print "[$TestDir][t$TID] Finished GUI";
# Compare with reference files (actually, completely fake at present, based on non-existent log file ...)
        $REFlogfn = &transpath("$REF/$logfn");
        print "[$TestDir][t$TID] Finished gui, comparing $logfn with $REFlogfn";
        system ("$DIFF_COMMAND $logfn $REFlogfn > $difffn");
        &dofldiff("Case_$TID","$logfn","$REF/$logfn",$fldopt);      # compare floating point numbers
        return;                 # go read next line from configuration file
        }

# ****************************************

sub domrggkktest {
        local($TID) = @_;
#
# purpose: exercise mrggkk
# arguments:
#               $TID = test identification that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# $fn1 = input to cut3d
# set default files names derived from $TID:
        $fn1 = "$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../telphon_3.in
        $infn = &transpath($fn1);
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = "t$TID.out";              # t80.out or t80.log
        $errfn = "t$TID.err";
        unlink($outfn,$errfn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting mrggkk";
        if ($OStype eq 'MacOS') {
                system ("$CODE_MRGGKK < $infn \xb7 $outfn ");
                }
        else {
		$REDIRECT_ERR = "2> $errfn";
                system ("$CODE_MRGGKK < $infn > $outfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
                }
        print "[$TestDir][t$TID] Finished mrggkk";
        return;                 # go read next line from configuration file
        }

# ****************************************
sub domrgtest {
	local($TID,$parms) = @_;
	local($fntodelete);
#
# purpose: exercise mrgddb
# arguments:
#	$TID = test identification that will be used as suffix to diff and prefix
# to log file names
# $parms = parameters string following the test ID; the format is:
#	[fn1=path] desc='short description' files=fn2,fn3[,...,fnLast] [opt=option]
#	path = output file name
#	fn2 ... fnLast = list of files to be mrgddb; variable $nfile will be
# set to the files count;
#	option is an option for fldiff, see setfnopt for the allowed values
#   fn1, description, $nfile and fn2 to fnLast will be put into the "files" file
# named "mrgddb.run"; this file will be used as input to mrgddb;
#   fn2 to fnLast can reside in the Test directory as normal or g-zipped files,
# they will be copied, and unzipped if necessary, to the working directory
# before the test and the copy will be deleted after
#   set default file name(s) derived from $TID:
	$fn1 = "t$TID.ddb.out";				# t14.ddb.out
	$descr = '';
# search description string first
	$ix = index ($parms,"desc='");
	$iy = $ix >= 0 ? index($parms,"'",$ix+6) : -1;
	if ($iy <= $ix) {
                print 'Error, description is missing or not surrounded by quotes';
                return;
                }
# check options left and right of description
	$descr = substr($parms,$ix+6,$iy - $ix -6);
	$subparms = substr($parms,0,$ix);
	@tokens = split(' ',$subparms);
	for ($i=0;$i <= $#tokens;$i++) {
		&setfnopt(@tokens[$i]);
		}
	$subparms = substr($parms,$iy+1);
	@tokens = split(' ',$subparms);
	for ($i=0;$i <= $#tokens;$i++) {
		&setfnopt(@tokens[$i]);
		}
	$nfile = $#mrgfiles +1;
	if ($debug >= 2) {
		print "FN1 $fn1 DESC $descr FILES $nfile @mrgfiles";
		print "fldflag $fldflag fldopt $fldopt";
		}
	if ($nfile <= 1) {
		print 'Error, files list is missing or not coma-separated';
		return;
		}
#	print "Case_$TID:";
# Nothing to do in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'mrgddb.run';		# files file for mrgddb
	$logfn = "t$TID.log";		# t14.log
	$errfn = "t$TID.err";		# t14.log
	$difffn = "diff.t$TID";		# diff.t14
	$outfn = &transpath($fn1);	# translate path name when required
	unlink($RUNfile,$difffn,$logfn,$errfn,$outfn);
# if a file fn2 ... fnLast does not exist, it is probably to be copied/unzipped from
# the parent directory
	@fntodelete = ();		# possible files to unzip/delete
	for ($n = 0; $n < $nfile; $n++) {
		if (! -f $mrgfiles[$n]) {
			@fntodelete = (@fntodelete, &copyunzip("$top_testdir/$TestDir/Input/$mrgfiles[$n]"));	# copy/unzip
			}
		}
	$fntodelete = &copyunzip("$top_testdir/$TestDir/Input/$fntounzip") if ($fntounzip ne '') ;	# copy/unzip if any
# Create "files" file with file names $fn1, ... one per line
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES $outfn;	# translate path name when required
	print FILES $descr;	# put description line
	print FILES $nfile;	# put files count
	for ($n = 0; $n < $nfile; $n++) {
		print FILES &transpath("$mrgfiles[$n]");	# translate path name when required
		}
	close(FILES);
	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] File $RUNfile created, starting mrgddb";
	if ($OStype eq 'MacOS') {
		system ("$CODE_MRGDDB < $RUNfile \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_MRGDDB < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
# Compare with reference files
	$REFoutfn = &transpath("$REF/$fn1");
	print "[$TestDir][t$TID] Finished mrgddb, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);	# compare floating point numbers
# check for work file to be deleted
	for ($n = 0; $n <= $#fntodelete; $n++) {
		system("$ERASE_COMMAND $fntodelete[$n]");	# delete file
		}
	return;			# go read next line from configuration file
	}

# ****************************************
sub domrgscrtest {
        local($TID) = @_;
#
# purpose: exercise mrgscr
# arguments:
#       $TID = test identification (1 to 3 characters string)
# internal file names:
#       $fn1 = input to mrgscr
#       $fn2 = output file name that will be compared with reference.
#       $fn3 will be used as suffix to diff and prefix to log file names
# set default files names derived from $TID:
        $fn1 = "$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../t79.in
        $fn2 = "t$TID.out";             # t79.out
        $fn3 = "t$TID";                 # t79
# set new file names or options according to the parameters
        for ($i=1;$i <= 9;$i++) {
                last if (@_[$i] eq '');
                &setfnopt(@_[$i]);
                }
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }

# Remove existing files related to this test case:
        $infn = &transpath($fn1);       # translate path name when required
        $outfn = &transpath($fn2);      # translate path name when required
        $errfn = $fn3.'.err';           # t79.err
        $difffn = 'diff.'.$fn3;         # diff.t79
        unlink($outfn,$logfn,$errfn,$difffn);

# Run the test case
        print "\n[$TestDir][t$TID] Starting mrgscr";
        if ($OStype eq 'MacOS') {
                system ("$CODE_MRGSCR < $infn \xb7 $logfn ");
                }
        else {
                $REDIRECT_ERR = "2> $errfn";
                system ("$CODE_MRGSCR < $infn > $outfn $REDIRECT_ERR");
                $REDIRECT_ERR = '2>&1';
                }

# Compare with reference files
        $REFoutfn = &transpath("$REF/$fn2");
        print "[$TestDir][t$TID] Finished mrgscr, comparing $outfn with $REFoutfn";
        system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
        &dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);  # compare floating point numbers
        return;                 # go read next line from configuration file
        }

# ****************************************
sub docut3dtest {
	local($TID,$p1,$p2,$p3,$p4,$p5) = @_;
#
# purpose: exercise cut3d
# arguments:
#		$TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
#		$p1,..$p5 = 1 to 5 parameters [re]defining some option or the file names (fn1 to fn4)
# the parameter format for a file name is: fnI=path
# the format for an option is: opt=-diflog (only allowed for now, see setfnopt)
# $fn1 = input to cut3d
# $fn2 = output file name that will be compared with reference.
# $fn3 = file name for the cut.in file
# $fn4 = file name for the xyz file
#   fn3 and fn4 can reside in the Test directory as normal or g-zipped files,
# they will be copied, and unzipped if necessary, to the working directory before
# the test and the copy will be deleted after
# set default files names derived from $TID:
	$fn1 = "$top_testdir/$TestDir/Input/t$TID.in";	# e.g. ../t77.in
	$fn2 = "t$TID.out";		# t77.out
	$fn3 = '';		# undefined
	$fn4 = '';		# undefined
# set new file names or option according to the parameters
	for ($i=1;$i <= 2;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#	print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
	$infn = &transpath($fn1);
	$difffn = "diff.t$TID";       # diff.t77
	$outfn = &transpath($fn2);      # t77.out
	$outdenfn = "t$TID.outden";# t77.outden
#	for cases v2 #80, v3 #59,61,63,64,66,67 set .log file to .out (not the usual log file)
	$logfn = $fldflag == $fldlog ? "t$TID.out" : "t$TID.log";		# t80.out or t80.log
	$errfn = "t$TID.err";
	unlink($logfn,$errfn,$difffn,$outfn,$outdenfn);
	if ($fn3 ne '') {
		$cutin = &copyunzip("$top_testdir/$TestDir/Input/$fn3");     # copy t78.cut.in
		$rc = system("$COPY_COMMAND $cutin cut.in");
		if ($rc != 0) {
			print "Error $rc copying file $cutin";
			return 32;
			}
		}
	$xyzin = &copyunzip("$top_testdir/$TestDir/Input/$fn4") if ($fn4 ne '') ;		# copy t78.xyz.in

	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] Starting cut3d";
	if ($OStype eq 'MacOS') {
		system ("$CODE_CUT3D < $infn \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_CUT3D < $infn > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
if ($fldflag eq $fldNorun){
print "[$TestDir][t$TID] No comparison made";
return;   #don't call fldiff
}

# Compare with reference files
	$REFoutfn = &transpath("$REF/$fn2");
	print "[$TestDir][t$TID] Finished cut3d, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);  # compare floating point numbers
	return;			# go read next line from configuration file
	}
# ****************************************
sub doaimtest {
	local($TID,$p1,$p2,$p3,$p4,$p5,$p6) = @_;
#
# purpose: exercise the atom-in-molecule Bader analysis (aim)
# arguments:
#		$TID = test identification (1 to 3 characters string)
#		$p1, ...,$p4 [,$p5] = 4 to 6 parameters [re]defining an option or the file
# names fn1 to fn5 that will be put into the "files" file named "aimfiles.run";
# this file will be used as input to aim
# $fn2 is the output file name that will be compared with reference
# $fn3 will be used as suffix to diff and prefix to log file names
# pseudo-potential file names can be specified for fn4, fn5 as psp1, psp2
# the parameters format is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
# set default files names derived from $TID:
	$fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";	# e.g. ../t57.in
	$fn2 = "t$TID".'i_DEN';		# t57i_DEN
	$fn3 = "t$TID";				# t57
	$fn4 = '';		# undefined
	$fn5 = '';		# undefined
	$psp[1] = '';		# undefined
	$psp[2] = '';		# undefined
# set new file names or options according to the parameters
	for ($i=1;$i <= 9;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
	$fn4 = $psp[1] if ($psp[1] ne '');
	$fn5 = $psp[2] if ($psp[2] ne '');
#	print "Case_$TID:";
# Test input file in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                system ("$CHKINABI -m aim $fn1 >> $CHKREPORT");
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'aimfiles.run';      # files file for aim
  $logfn = $fn3.'.log';           # t57.log
  $errfn = $fn3.'.err';           # t57.err
  $difffn = 'diff.'.$fn3;         # diff.t57
  $outfn = $fn3.'.out';           # t57.out
  unlink($RUNfile,$difffn,$logfn,$errfn,$outfn);
  $fnglob = "$fn3.$MULTI";            # t57.*
  system("$ERASE_COMMAND $fnglob");		# unlink < t57.* >
# Create "files" file with file names $fn1, ... one per line
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES &transpath($fn1);   # translate path name when required
	print FILES $fn2;
	print FILES $fn3;
	print FILES &transpath($fn4);   # translate path name when required
	print FILES &transpath($fn5) if ($fn5 ne '');   # add 5th filename if any
	close(FILES);
	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] File $RUNfile created, starting aim";
	if ($OStype eq 'MacOS') {
		system ("$CODE_AIM < $RUNfile \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_AIM < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
# Compare with reference files
	$REFoutfn = &transpath("$REF/$outfn");
	print "[$TestDir][t$TID] Finished aim, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$outfn","$REF/$outfn",$fldopt);  # compare floating point numbers
	return;              # go read next line from configuration file
	}
# ****************************************
sub doconductitest {
        local($TID,$p1,$p2,$p3) = @_;
#
# purpose: exercise conducti
# arguments:
#		$TID = test identification (1 to 3 characters string)
#		$p1,$p2,$p3 = parameters [re]defining the file names fn1 to fn3 that will
# be put into the "files" file named "conducti.run";
# the parameters format is: fnI=path
#		$fn1 = input to conducti
#		$fn2 = output file name that will be compared with reference.
# 	$fn3 will be used as suffix to diff and prefix to log file names
# set default files names derived from $TID:
	$fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";	# e.g. ../t79.in
	$fn2 = "t$TID.out";		# t79.out
	$fn3 = "t$TID";			# t79
# set new file names or options according to the parameters
	for ($i=1;$i <= 9;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#	print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'conducti.run'; # files file for conducti
	$logfn = $fn3.'.log';           # t79.log
	$errfn = $fn3.'.err';           # t79.err
	$difffn = 'diff.'.$fn3;         # diff.t79
	$outfn = &transpath($fn2);      # translate path name when required
	unlink($RUNfile,$difffn,$logfn,$errfn,$outfn);

# Create "files" file with file names $fn1,2
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES &transpath($fn1);   # translate path name when required
	print FILES &transpath($fn3);   # translate path name when required
	close(FILES);

	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] Starting conducti";
	if ($OStype eq 'MacOS') {
		system ("$CODE_CONDUCTI < $RUNfile \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_CONDUCTI < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}

# Compare with reference files
	$REFoutfn = &transpath("$REF/$fn2");
	print "[$TestDir][t$TID] Finished conducti, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);  # compare floating point numbers
	return;                 # go read next line from configuration file
	}

# ****************************************
sub dofftproftest {
        local($TID) = @_;
#
# purpose: use fftprof to test FFT routines at the unitary level
# arguments:
#       $TID = test identification (1 to 3 characters string)
# internal file names:
#       $fn1 = input to fftprof
#       $fn2 = output file name that will be compared with reference.
#       $fn3 will be used as suffix to diff and prefix to log file names
# set default files names derived from $TID:
        $fn1 = "$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../t79.in
        $fn2 = "t$TID.out";             # t79.out
        $fn3 = "t$TID";                 # t79
# set new file names or options according to the parameters
        for ($i=1;$i <= 9;$i++) {
                last if (@_[$i] eq '');
                &setfnopt(@_[$i]);
                }
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }

# Remove existing files related to this test case:
        $infn = &transpath($fn1);       # translate path name when required
        $outfn = &transpath($fn2);      # translate path name when required
        $errfn = $fn3.'.err';           # t79.err
        $difffn = 'diff.'.$fn3;         # diff.t79
        unlink($outfn,$logfn,$errfn,$difffn);

# Run the test case
        print "\n[$TestDir][t$TID] Starting fftprof";
        if ($OStype eq 'MacOS') {
                system ("$CODE_FFTPROF < $infn \xb7 $logfn ");
                }
        else {
                $REDIRECT_ERR = "2> $errfn";
                system ("$CODE_FFTPROF < $infn > $outfn $REDIRECT_ERR");
                $REDIRECT_ERR = '2>&1';
                }

# Compare with reference files
        $REFoutfn = &transpath("$REF/$fn2");
        print "[$TestDir][t$TID] Finished fftprof, comparing $outfn with $REFoutfn";
        system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
        &dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);  # compare floating point numbers
        return;                 # go read next line from configuration file
        }

# ****************************************
sub doujdettest {
        local($TID) = @_;
#
# purpose: exercise ujdet
# arguments:
#               $TID = test identification (1 to 3 characters string)
#               $p1,$p2,$p3 = parameters [re]defining the file names fn1 to fn3 that will
# be put into the "files" file named "ujdet.run";
# the parameters format is: fnI=path
#               $fn1 = input to ujdet
#               $fn2 = output file name that will be compared with reference.
#       $fn3 will be used as suffix to diff and prefix to log file names
# set default files names derived from $TID:
        $fn2 = "ujdet.out";
        $fn3 = "t$TID";
        $fn4 = "t$TID.out";
# set new file names or options according to the parameters
        for ($i=1;$i <= 9;$i++) {
                last if (@_[$i] eq '');
                &setfnopt(@_[$i]);
                }
#       print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $RUNfile = 'ujdet.run'; # files file for conducti
        $logfn = $fn3.'.log';           # t79.log
        $errfn = $fn3.'.err';           # t79.err
        $difffn = 'diff.'.$fn3;         # diff.t79
        $outfn = &transpath($fn2);      # translate path name when required
        unlink($RUNfile,$difffn,$logfn,$errfn,$outfn);

# Run the test case
        print "\n[$TestDir][t$TID] Starting ujdet";
        if ($OStype eq 'MacOS') {
                system ("$CODE_UJDET \xb7 $logfn ");
                }
        else {
                $REDIRECT_ERR = "2> $errfn";
                system ("$CODE_UJDET > $logfn $REDIRECT_ERR");
                $REDIRECT_ERR = '2>&1';
                }

# Compare with reference files
        $REFoutfn = &transpath("$REF/$fn4");
        print "[$TestDir][t$TID] Finished ujdet, comparing $outfn with $REFoutfn";
        system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
        &dofldiff("Case_$TID","$fn2","$REF/$fn4",$fldopt);  # compare floating point numbers
        return;                 # go read next line from configuration file
        }

# ****************************************
sub dooptictest {
        local($TID,$p1,$p2,$p3) = @_;
#
# purpose: exercise optic
# arguments:
#               $TID = test identification (1 to 3 characters string)
#               $p1,$p2,$p3 = parameters [re]defining the file names fn1 to fn3 that will
# be put into the "files" file named "optic.run";
# the parameters format is: fnI=path
#               $fn1 = input to optic
#               $fn2 = output file name that will be compared with reference.
#       $fn3 will be used as suffix to diff and prefix to log file names, and also for roots for temporaries
# set default files names derived from $TID:
        $fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";   # e.g. ../t79.in
        $fn2 = "t$TID.out";             # t79.out
        $fn3 = "t$TID";                         # t79
# set new file names or options according to the parameters
        for ($i=1;$i <= 9;$i++) {
                last if (@_[$i] eq '');
                &setfnopt(@_[$i]);
                }
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $RUNfile = 'optic.run'; # files file for optic
        $logfn = $fn3.'.log';           # t79.log
        $errfn = $fn3.'.err';           # t79.log
        $difffn = 'diff.'.$fn3;         # diff.t79
        $fn3linopt = $fn3.'_0001_0001-linopt.out';
        $fn3chitotre = $fn3.'_0001_0002_0003-ChiTotRe.out';
        $fn3chitotim = $fn3.'_0001_0002_0003-ChiTotIm.out';
        $outfnl = &transpath($fn3linopt);        # translate path name when required
        $outfnr = &transpath($fn3chitotre);      # translate path name when required
        $outfni = &transpath($fn3chitotim);      # translate path name when required
        unlink($RUNfile,$difffn,$logfn,$errfn,$outfnl,$outfnr,$outfni);

# Create "files" file with file names $fn1,2,3
        open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
        print FILES &transpath($fn1);   # translate path name when required
        print FILES &transpath($fn2);   # translate path name when required
        print FILES &transpath($fn3);   # translate path name when required
        close(FILES);

        return if ($debug >= 3);
# Run the test case
        print "\n[$TestDir][t$TID] Starting optic";
        if ($OStype eq 'MacOS') {
                system ("$CODE_OPTIC < $RUNfile \xb7 $logfn ");
                }
        else {
		$REDIRECT_ERR = "2> $errfn";
                system ("$CODE_OPTIC < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
                }

# Compare with reference files
        $REFoutfn = &transpath("$REF/$fn3linopt");
        print "[$TestDir][t$TID] Finished optic, comparing $outfnl with $REFoutfn";
        system ("$DIFF_COMMAND $outfnl $REFoutfn > $difffn");
        &dofldiff("Case_$TID : linopt ","$fn3linopt","$REF/$fn3linopt",$fldopt);  # compare floating point numbers

        $REFoutfn = &transpath("$REF/$fn3chitotre");
        if(-f $REFoutfn) {
          print "[$TestDir][t$TID] Finished optic, comparing $outfnr with $REFoutfn";
          system ("$DIFF_COMMAND $outfnr $REFoutfn > $difffn");
          &dofldiff("Case_$TID : chitotre ","$fn3chitotre","$REF/$fn3chitotre",$fldopt);  # compare floating point numbers
                         }

        $REFoutfn = &transpath("$REF/$fn3chitotim");
        if(-f $REFoutfn) {
          print "[$TestDir][t$TID] Finished optic, comparing $outfni with $REFoutfn";
          system ("$DIFF_COMMAND $outfni $REFoutfn > $difffn");
          &dofldiff("Case_$TID : chitotim ","$fn3chitotim","$REF/$fn3chitotim",$fldopt);  # compare floating point numbers
                         }

        return;                 # go read next line from configuration file
        }

# ****************************************
sub dolwftest {
	local($TID,$p1,$p2,$p3,$p4,$p5) = @_;
#
# purpose: exercise lwf
# arguments:
#		$TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
#		$p1,..$p5 = 1 to 5 parameters [re]defining some option or the file names (fn1 to fn4) that
# will be put into the "files" file named "lwf.run";
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#	$fn1 = main input to lwf
# $fn2 = moldyn input to lwf
# $fn3 = output file name that will be compared with reference.
# $fn4 = auxiliary output file name
# set default files names derived from $TID:
	$fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";		# e.g. ../t88.in
	$fn2 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.moldyn.in";	# ../t88.moldyn.in
	$fn3 = "t$TID.out";			# t88.out
	$fn4 = "t$TID.wandata";		# t88.wandata
# set new file names or options according to the parameters
	for ($i=1;$i <= 5;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#	print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'lwf.run';           # files file for lwf
	$infn = &transpath($fn1);
	$logfn = "t$TID.log";		# t88.log
	$errfn = "t$TID.err";		# t88.err
	$difffn = 'diff.t'.$TID;       # diff.t88
	$outfn = &transpath($fn3);      # t88.out
	$wandatafn = &transpath($fn4);      # t88.wandata
	unlink($RUNfile,$logfn,$errfn,$difffn,$outfn,$wandatafn);
# Create "files" file with file names $fn1, ... one per line
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES &transpath($fn1);   # translate path name when required
	print FILES &transpath($fn2);
	print FILES $fn3;
	print FILES $fn4;
	close(FILES);
	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] Starting lwf";
	if ($OStype eq 'MacOS') {
		system ("$CODE_LWF < $RUNfile \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_LWF < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}

# Compare with reference files
	$REFoutfn = &transpath("$REF/$fn3");
	print "[$TestDir][t$TID] Finished lwf, comparing $outfn with $REFoutfn";
	system ("$DIFF_COMMAND $outfn $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn3","$REF/$fn3",$fldopt);  # compare floating point numbers
	return;                 # go read next line from configuration file
	}
# ****************************************
sub dodocchk {
        local($TID) = @_;
#
# purpose: check the documentation
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
        $fn2 = "t$TID.log";                 # t51.log
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $logfn = &transpath($fn2);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting docchk";
        if ($OStype eq 'MacOS') {
                system ("$CODE_DOCCHK \xb7 $logfn ");
                }
        else {
		print "$CODE_DOCCHK > $logfn $REDIRECT_ERR" ;
                system ("$CODE_DOCCHK > $logfn $REDIRECT_ERR");
                }
        system ("grep -v SUCCESS $logfn > $outfn");
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished docchk, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub docheck_forbidden {
        local($TID) = @_;
#
# purpose: check the existence of forbidden access to standard output
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting check_forbidden";
        if ($OStype eq 'MacOS') {
                system ("$CODE_CHECK_UNALLOWED \xb7 $outfn ");
                }
        else {
                print "$CODE_CHECK_UNALLOWED > $outfn $REDIRECT_ERR" ;
                system ("$CODE_CHECK_UNALLOWED > $outfn $REDIRECT_ERR");
                }
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished check_forbidden, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub docheck_forbidden_in_doc_yml {
        local($TID) = @_;
#
# purpose: check the existence of forbidden strunk in doc yml files
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting check_forbidden_in_doc_yml";
        if ($OStype eq 'MacOS') {
                system ("$CODE_CHECK_FORBIDDEN_IN_DOC_YML \xb7 $outfn ");
                }
        else {
                print "$CODE_CHECK_FORBIDDEN_IN_DOC_YML > $outfn $REDIRECT_ERR" ;
                system ("$CODE_CHECK_FORBIDDEN_IN_DOC_YML > $outfn $REDIRECT_ERR");
                }
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished check_forbidden_in_doc_yml, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub docheck_ascii {
        local($TID) = @_;
#
# purpose: check the existence of non-ASCII characters in source files
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting check_ascii";
        if ($OStype eq 'MacOS') {
                system ("$CODE_CHECK_ASCII \xb7 $outfn ");
                }
        else {
                print "$CODE_CHECK_ASCII > $outfn $REDIRECT_ERR" ;
                system ("$CODE_CHECK_ASCII > $outfn $REDIRECT_ERR");
                }
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished check_ascii, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub docheck_inlined_macros {
        local($TID) = @_;
#
# purpose: check for the presence of inlined macros in source files
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting check_inlined_macros";
        if ($OStype eq 'MacOS') {
                system ("$CODE_CHECK_INLINED_MACROS \xb7 $outfn ");
                }
        else {
                print "$CODE_CHECK_INLINED_MACROS > $outfn $REDIRECT_ERR" ;
                system ("$CODE_CHECK_INLINED_MACROS > $outfn $REDIRECT_ERR");
                }
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished check_inlined_macros, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub dochkwarnings {
        local($TID,$p1) = @_;
#
# purpose: check the warnings
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
        $fn2 = "t$TID.log";                 # t51.log
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $logfn = &transpath($fn2);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting WARNCHK";
        if ($OStype eq 'MacOS') {
                system ("$CODE_WARNCHK $p1 \xb7 $logfn ");
                }
        else {
		print "$CODE_WARNCHK $p1 > $logfn $REDIRECT_ERR" ;
                system ("$CODE_WARNCHK $p1 > $logfn $REDIRECT_ERR");
                if ( $p1 eq '3' or $p1 eq '4') {
                    print "$SORTED_CMD $logfn" ;
                    system ("$SORTED_CMD $logfn")
                    }
                }
        system ("grep -v SUCCESS $logfn > $outfn");
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished WARNCHK, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub dochkinpvars {
        local($TID) = @_;
#
# purpose: check the documentation
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";                 # t51.out
        $fn2 = "t$TID.log";                 # t51.log
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath($fn1);
        $logfn = &transpath($fn2);
        $difffn = "diff.t$TID";      # diff.t51
        unlink($RUNfile,$outfn,$logfn,$difffn);
# Run the test case
        print "\n[$TestDir][t$TID] Starting dochkinpvars";
        if ($OStype eq 'MacOS') {
                system ("$CODE_CHKINPVARS \xb7 $logfn ");
                }
        else {
                print "$CODE_CHKINPVARS > $logfn $REDIRECT_ERR" ;
                system ("$CODE_CHKINPVARS > $logfn $REDIRECT_ERR");
                }
        system ("grep -v SUCCESS $logfn > $outfn");
# Compare with reference files
        $REFfn = &transpath("$REF/$fn1");
        print "[$TestDir][t$TID] Finished docheckinputvars, comparing $outfn with $REFfn";
        system ("$DIFF_COMMAND $outfn $REFfn > $difffn");
        &dofldiff("Case_$TID","$fn1","$REF/$fn1",$fldopt);  # compare
        return;                 # go read next line from configuration file
        }

# ****************************************
sub doband2epstest {
	local($TID,$p1,$p2,$p3,$p4,$p5) = @_;
#
# purpose: exercise band2eps
# arguments:
#		$TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
#		$p1,..$p5 = 1 to 5 parameters [re]defining some option or the file names (fn1 to fn4) that
# will be put into the "files" file named "band2eps.run";
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#	$fn1 = input to band2eps
#	$fn2 = output file (eps) name that will be compared with reference.
#	$fn3 = _freq file
#	$fn4 = _displ file
# set default files names derived from $TID:
	$fn1 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in";		# e.g. ../t51.in
	$fn2 = "t$TID.out.eps";			# t51.out.eps
	$fn3 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in_freq";	# ../t51.in_freq
	$fn4 = "$CYGWIN$top_testdir/$TestDir/Input/t$TID.in_displ";	# ../t51.in_displ
# set new file names or options according to the parameters
	for ($i=1;$i <= 5;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#	print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'band2eps.run';     # files file for band2eps
	$infn = &transpath($fn1);
	$logfn = "t$TID.log";     # t51.log
	$errfn = "t$TID.err";     # t51.log
	$difffn = "diff.t$TID";      # diff.t51
	$epsfn = &transpath($fn2);
	$freqfn = &transpath($fn3);    # .freq (input)
	$displfn = &transpath($fn4);   # .displ (output)
	unlink($RUNfile,$logfn,$errfn,$difffn,$epsfn);
# Create "files" file with file names $fn1, ... one per line
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	print FILES $infn;   # translate path name when required
	print FILES $epsfn;
	print FILES $freqfn;
	print FILES $displfn;
	close(FILES);
	return if ($debug >= 3);
# Run the test case
	print "\n[$TestDir][t$TID] File $RUNfile created, starting band2eps";
	if ($OStype eq 'MacOS') {
		system ("$CODE_BAND2EPS < $RUNfile \xb7 $logfn ");
		}
	else {
		$REDIRECT_ERR = "2> $errfn";
		system ("$CODE_BAND2EPS < $RUNfile > $logfn $REDIRECT_ERR");
    		$REDIRECT_ERR = '2>&1';
		}
# Compare with reference files
	$REFepsfn = &transpath("$REF/$fn2");
	print "[$TestDir][t$TID] Finished band2eps, comparing $epsfn with $REFepsfn";
	system ("$DIFF_COMMAND $epsfn $REFepsfn > $difffn");
	&dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);  # compare floating point numbers
	return;                 # go read next line from configuration file
	}
# ****************************************
sub domacroavetest {
	local($TID,$p1,$p2,$p3,$p4,$p5) = @_;
#
# purpose: exercise macroave
# arguments:
#		$TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
#		$p1,..$p5 = 1 to 5 parameters [re]defining some option or the file names (fn1 to fn4) that
# will be put into the "files" file named "macroave.in";
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#	$fn1 = input to macroave
#	$fn2 = first output file name that will be compared with reference.
#	$fn3 = second output file name that will be compared with reference.
#	$fn4 is also test number that will be used as suffix to diff
# set default files names derived from $TID:
	$fn1 = "$top_testdir/$TestDir/Input/t$TID.in";		# e.g. ../t41.in
	$fn2 = "t$TID".'_DEN.MAV';		# t41_DEN.MAV
	$fn3 = "t$TID".'_DEN.PAV';		# t41_DEN.PAV
	$fn4 = "t$TID";			# t41
# set new file names or options according to the parameters
	for ($i=1;$i <= 5;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}

#	print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
	$RUNfile = 'macroave.in'; # files file for macroave
	$difffn = 'diff.'.$fn4;         # diff.t41
	$outfn1= &transpath($fn2);      # translate path name when required
	$outfn2= &transpath($fn3);      # translate path name when required
	unlink($RUNfile,$difffn,$outfn1,$outfn2);

# Create "macroave.in" file from file $fn1
	$rc1 = system("$COPY_COMMAND $fn1 $RUNfile");
	print "\n[$TestDir][t$TID] Case_$TID: file $fn1 copied to $RUNfile" if ($rc1 == 0);

	return if ($debug >= 3);
# Run the test case
	print "[$TestDir][t$TID] Starting macroave";
	system ("$CODE_MACROAVE");

# Compare with reference files
	$REFoutfn = &transpath("$REF/$fn2");
	print "[$TestDir][t$TID] Finished macroave, comparing $outfn1 with $REFoutfn";
	system ("$DIFF_COMMAND $outfn1 $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn2","$REF/$fn2",$fldopt);  # compare floating point numbers
	$REFoutfn = &transpath("$REF/$fn3");
	print "[$TestDir][t$TID] Finished macroave, comparing $outfn2 with $REFoutfn";
	system ("$DIFF_COMMAND $outfn1 $REFoutfn > $difffn");
	&dofldiff("Case_$TID","$fn3","$REF/$fn3",$fldopt);  # compare floating point numbers
	return;                 # go read next line from configuration file
	}
# ****************************************
sub doserieX {
	local($TID,$p1,$p2,$p3,$p4,$p5) = @_;
#
# purpose: run cpu tests of the sequential version of abinit
#
# arguments:
#		$TID = test identification (1 letter X and 1 or 2 other characters) that will be used as suffix to several
# file names
#		$p1,..$p5 = 1 to 5 parameters [re]defining some option or the file names (fn1,fn2)
# the parameter format for a file name is: fnI=path
# the parameter format for the FFT algorithm is fftalg=abc
# the parameter format for the other options is: opt=option ; see setfnopt
#	$fn1 = name of common input file for the serie
# $fn2 = name of the reference file for the serie
# set default files names derived from $TID:
	$X = substr($TID,0,1);		# get serie letter [A-D] from the test id
	$fn1 = "input$X";			# e.g. inputA
	$fn2 = "report$X".'3';			# reportA3
	$FFTalg = '';		# undefined
# set new file names or options according to the parameters
	for ($i=1;$i <= 5;$i++) {
		last if (@_[$i] eq '');
		&setfnopt(@_[$i]);
		}
#	print "Case_$TID:";
# Build input file
	$infname = &transpath("$top_testdir/$TestDir/Input/$fn1");	# translate path if necessary
	$inputXn = "input$TID";	# file name of the input file for this serie e.g. inputA1
	$rc = system("$COPY_COMMAND $infname $inputXn");	# copy common part of file
  if ($rc != 0) {
  	print "Error $rc copying file $infname";
		return 32;
		}
# check for options
	if ($FFTalg ne '' || $option != 0) {
		open (INFILE,">>$inputXn") || die "Unable to open $inputXn";
		print INFILE " fftalg $FFTalg" if ($FFTalg ne '');		# append FFT algorithm if necessary
		print INFILE ' istwfk 1' if ($option == $optistwfk1);	 # append istwfk values if it applies
		print INFILE ' istwfk 1 1' if ($option == $optistwfk11);	 # append istwfk values if it applies
		close INFILE;
		}
# Establish test subdirectory
	$SERIE_DIR = $TID;
	if (! -e $SERIE_DIR || ! -d $SERIE_DIR) {
		mkdir ($SERIE_DIR,0755);		# Mode 0755 ignored under DOS-Windows
		print "Subdirectory $SERIE_DIR created";
		}
	else {
		print "Cannot create subdirectory, $SERIE_DIR already exists";
		}
# Change to the test subdirectory
	print "cd $SERIE_DIR";
	chdir ("$SERIE_DIR");
  $WKdir = &transpath("$WORK_DIR/$SERIE_DIR");
	print "Working directory changed : now $WKdir \n";
	$RUNfile = 'filnames.run';	# files file for abinit
	$reportfn = 'report';		# default report file name for getrep
	unlink ($RUNfile,$reportfn);	# erase old files
	system ("$ERASE_RECURSIVE run*");	# erase recursively run* directories

# copy "files" file from upper directory and translate file pathes
	$ipath = &transpath('../../files');
	open(FILEI,"<$ipath") || die "Unable to open file files for test $TID";
	open(FILES,">$RUNfile") || die "Unable to open FILES for test $TID";
	$inputfn = <FILEI>;	# first line from files is also input file for abinit
	chop $inputfn;		# drop NL
	$inputfn = &transpath($inputfn);
	print FILES $inputfn;	# rewrite translated path of input
	while (<FILEI>) {	# read until end of file
		chop $_;	# drop NL
		print FILES &transpath($_);	# rewrite other translated pathes
		}
	close FILEI;
	close FILES;
# copy input file
	unlink ($inputfn);		# erase old input file
	$inputpath = &transpath("../$inputXn");	# translate path if necessary
	$rc = system("$COPY_COMMAND $inputpath $inputfn");
	if ($rc != 0) {
		print "Error $rc copying file $inputXn";
		return 32;
		}
	return if ($debug >= 3);
	print "Running serie";
	$GETREP = $PERL.&perlpath(" $top_testdir/scripts/make-cpu-report.pl");	# command to call getrep perl script
	$REFDIR = &transpath("$top_testdir/$TestDir/Refs");	# reference directory
# DEFAULT NAMES for the files, input, log, ... files for abinit and getrep
	$FNinput = 'input';
	$LOGFILE = 'log';
	$FNout= 'output';	# abinit output file
	$DIRPRFX = 'run';	# prefix for the tests subdirectories
	$REPlog = 'getrep.log';	# log file name for getrep
# Since DOS does not support the metacharacter * in the directory part of a fileid,
# the following variable will hold the list of files to be processed by the getrep script:
$REPfiles = '';		 # empty list for now

# A number of keywords, and the different values must be specified as input
# for the different runs.
# DEFAULT NAMES AND VALUES of the keywords for series A-D

	if ($X eq 'A') {		# This particular version is for the small FFT tests
		@keywd=('ecut','ngfft');	# 2 keywords, 6 runs
		@data = ('4.9','3*20',
			 '7'  ,'3*24',
			'10.8','3*30',
			'12.5','3*32',
			'15.9','3*36',
			'19.5','3*40');
		}

	elsif ($X eq 'B') {	# This particular version is for the big FFT tests
		@keywd=('ecut','ngfft');	# 2 keywords, 7 runs
		@data = ('28' ,'3*48',
			 '44' ,'3*60',
			 '50' ,'3*64',
			 '63' ,'3*72',
			 '78' ,'3*80',
			 '98' ,'3*90',
			'112' ,'3*96');
		}

	elsif ($X eq 'C') {	# This particular version is for the scaling tests
# 6 keywords, 4 runs
		@keywd=('acell',          'kpt',  'natom','nband','typat',	'xred');
		@data=( '3*10.00',     '.25 .25 .25' ,'2' ,'4'	 ,'2*1',	'5*.0 .4',
			'2*10.0 20.0', '.25 .25 .5'  ,'4' ,'8'	 ,'4*1', '5*.0 .2 2*.0 .5 2*.0 .7',
			'10.0 2*20.0', '.25 .5 .5'   ,'8','16'	 ,'8*1', '5*.0 .2 2*.0 .5 2*.0 .7  .0 .5 .0  .0 .5 .2 .0 2*.5 .0 .5 .7',
			'3*20.0',	'.5 .5 .5',  '16','32'	 ,'16*1','5*.0 .2 2*.0 .5 2*.0 .7  .0 .5 .0  .0 .5 .2 .0 2*.5 .0 .5 .7 .5 2*.0 .5 .0 .2 .5 .0 .5 .5 .0 .7  2*.5 .0  2*.5 .2 3*.5 2*.5 .7' );
		}

	elsif ($X eq 'D') {	# This particular version is for the FFT tests with
# unusual values
		@keywd=('ecut','ngfft');	# 2 keywords, 7 runs
		@data =( '7.1','3*25',
			 '8.2','3*27',
			'23'  ,'3*45',
			'30'  ,'3*50',
			'35'  ,'3*54',
			'67'  ,'3*75',
			'78.5','3*81' );
		}
	else {
		print "Unrecognized serie identifier $X";
		return;
		}
# Here begins the code : IT IS VERY HAZARDOUS TO MODIFY IT !
# ***********************************************************
	print "Here begins the execution of serie $X tests, $number_kw keywords, $number_run runs";

	$number_kw = $#keywd + 1;	# keywords count
	$number_run = ($#data + 1) / $number_kw;	# runs count
	$ix_data = 0;			# initialize index into $data array
	for ($ix_run = 1; $ix_run <= $number_run; $ix_run ++) {		# Loop on all tests
# Change the working directory
		$DIR = $DIRPRFX.$ix_run;
		$rc = system("$ERASE_RECURSIVE $DIR");
		mkdir ($DIR,0755);		# Mode 0755 ignored under DOS-Windows
		$cmd = "$COPY_COMMAND $RUNfile $DIR";
		$rc = system($cmd);
		if ($rc != 0) {
			print "Error $rc copying $RUNfile";
			return 100;
			}
		print "File $RUNfile copied";
		print "cd $DIR";
		chdir ("$DIR");
		$WKdir = &transpath("$WORK_DIR/$SERIE_DIR/$DIR");
		print "Working directory changed : now $WKdir \n";
# Read the file $FNinput in upper directory and build new file in current one;
# new values will be substituted in lines beginning by keywords from $keywd
		$upper_input = &transpath("../$FNinput");	# translate path if necessary
		open (INFILE,"<$upper_input") || die "Unable to open $upper_input";
		open (OUT,">$FNinput");		# open file named "input" in write mode
		while (<INFILE>) {
			chop $_;		# drop NL
			($inkeyword,$word1,$remain) = split(' ',$_);
			for ($ix_kwd = 0; $ix_kwd < $number_kw; $ix_kwd ++) {
				if ($inkeyword eq $keywd[$ix_kwd]) {
# point to corresponding new value in data array and build new line:
					$ix_value = $ix_data+$ix_kwd;
					$line = "$keywd[$ix_kwd] $data[$ix_value]";
					print OUT $line;	# write updated line to file OUT
					print "File $FNinput changed for: $line";
					last;
					}
				}	# End of the for loop on keywords
# if the keyword from $upper_input file has not been found in the $keywd array, the
# previous for loop will end with an highest index value; write input line unchanged
# to file OUT:
			print OUT $_ if ($ix_kwd >= $number_kw);
			}	# end of the while loop to read INFILE
		close INFILE;
		close OUT;
# Bump index to new values in the data array for the next run:
		$ix_data += $number_kw;
		print "Keywords processing is finished";
# Execute one test
		print "Calling the abinit code, sequential version";
		$cmd = "$CODE_CPU < $RUNfile > $LOGFILE";
		print $cmd;
		$rc  = system($cmd);
		$REPfiles .= &perlpath("$DIR/$FNout ");	# add output file to getrep list
		print "Done with $DIR";
		system("$ERASE_COMMAND *_WFK");
# Change the working directory to parent
		print 'cd ..';
	 	chdir ('..');
		print "Working directory changed\n";
		}		# End of the loop on the tests

# Now, execute the getrep script
	print "Running getrep";
	$cmd = "$GETREP $REPfiles > $REPlog $REDIRECT_ERR";
	print $cmd;
	$rc = system($cmd);

# Compare the result with the reference file
	print 'Compare results for energy with previous results';
	unlink("diff_$SERIE_DIR");
	$REFFILE = &transpath("$REFDIR/$fn2");	# reference file path
	$cmd = "$DIFF_COMMAND $reportfn $REFFILE > diff_$SERIE_DIR";
	print $cmd;
	system($cmd);

# Clean the directory
	system("$ERASE_COMMAND cpu*");
# Change the working directory to parent
	print 'cd ..';
	chdir ('..');
	print "Working directory changed";
# Delete input file for the current serie here
	unlink ("$inputXn");
	return;			# go read next line from configuration file
	}
# ****************************************
sub docopy {
	local($TID,$p1,$p2) = @_;
#
# purpose: copy file(s)
# arguments:
#		$TID = test identification (1 to 3 characters string)
#		$p1 = name of source file(s) to be copied
#		$P2 = destination file or directory
# Remove the $CYGWIN prefix for system commands like copy diff...
        if ( $p1 =~ m/\/cygwin(.*)/ ) { $p1 = $1; };
#
	$p1T = &transpath($p1);	# translate path name when required
	$cmd = "$COPY_COMMAND $p1T $p2";
	print "COPY: $cmd" if ($debug >= 2);
	$rc = system($cmd);
	if ($rc == 0) {
		print "\n[$TestDir][t$TID] File $p1 copied to $p2";
		}
	else {
		print "\n[$TestDir][t$TID] Error $rc copying $p1";
		}
	return;			# go read next line from configuration file
	}
# ****************************************
sub doerase {
	local($TID,$p1) = @_;
#
# purpose: erase file(s)
# arguments:
#		$TID = test identification (1 to 3 characters string)
#		$p1 = name of file(s) to be erased
	$cmd = "$ERASE_COMMAND $p1";
	print "ERASE: $cmd" if ($debug >= 2);
	$rc = system($cmd);
	if ($rc == 0) {
		print "\n[$TestDir][t$TID] File $p1 erased";
		}
	else {
		print "\n[$TestDir][t$TID] Error $rc erasing $p1";
		}
	return;			# go read next line from configuration file
	}
# ****************************************
sub doncd {
        local($TID,$p1,$p2) = @_;
#
# purpose : use the ncdump program to generate ASCII file from a netcdf file
# arguments :
#               $TID = test identification (1 to 3 characters string)
#               $p1 = netcdf file
#               $p2 = ASCII file
        system("$CODE_NCDUMP $p1 > $p2");
        print "creating ASCII file";

# Compare with reference files
        $REFoutfn = &transpath("$REF/$p2");
        print "Comparing $outfn with $REFoutfn";
        system ("$DIFF_COMMAND $p2 $REFoutfn > $difffn");
        &dofldiff("Case_$TID","$p2","$REF/$p2",$fldopt);  # compare floating point numbers
        return;                 # go read next line from configuration file
        }
# ****************************************
sub doreport {
        local($TID) = @_;
#
# purpose: exercise lwf
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
#               $p1,..$p5 = 1 to 5 parameters [re]defining some option or the file names (fn1 to fn4) that
# will be put into the "files" file named "lwf.run";
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = main input to lwf
# $fn2 = moldyn input to lwf
# $fn3 = output file name that will be compared with reference.
# $fn4 = auxiliary output file name
# set default files names derived from $TID:
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
	if ( $TestDir eq "buildsys" ) {
		return;
		}

# Prepare filenames
        $fn1 = "$top_testdir/$TestDir/Input/report.in";
        $reportref = &transpath($fn1);
        $fn2 = "fldiff.report" ;
        $fldiffreport = &transpath($fn2);
        $fn3 = "summaryfldiff" ;
        $summaryfldiff = &transpath($fn3);
        $fn4 = "report" ;
        $shortreport = &transpath($fn4);
# Create summary of tests
        unlink($summaryfldiff);
#        system ("grep 'Summary' $fldiffreport > $summaryfldiff");
        system ("grep --no-filename 'Summary' $fldiffreport > $summaryfldiff");
# Compare with reference files
        print "Prepare the short report";
        unlink($shortreport);
        &reportdiff("$summaryfldiff","$reportref","$shortreport");  # compare floating point numbers
        return;                 # go read next line from configuration file
	}

# ****************************************
sub dochkstatus {
        local($TID,$p1) = @_;
	if ( $debug >= 2) { print "Test $TID $p1"; }
#
# purpose: check the warnings
# arguments:
#               $TID = test identification (1 to 3 characters string) that will be used as suffix to diff
# and prefix to log file names
# the parameter format for a file name is: fnI=path
# the parameter format for an option is: opt=option ; see setfnopt
#       $fn1 = output file that will be compared with reference.
#       $fn2 = log
# set default files names derived from $TID:
        $fn1 = "t$TID.out";
        $fn2 = "t$TID.err";
        $fn3 = "report";
#        print "Case_$TID:";
# Nothing done in case 'MACH' is set to chkinabi
        if ($MACH eq 'chkinabi') {
                return ;
                }
# Remove existing files related to this test case:
        $outfn = &transpath("$abinit_outdir/$TestDir/$WORK_DIR/$fn1");
        $logfn = &transpath("$abinit_outdir/$TestDir/$WORK_DIR/$fn2");
        unlink($outfn,$logfn);
        $CODE_STATCHK = &transpath("$p1");
        $CODE_STATCHK = "$timeout_cmd$CODE_STATCHK";
# Run the test case
        print "\n[$TestDir][t$TID] Starting STATCHK";
        $my_path = $ENV{"PWD"};
        chdir($abinit_srcdir);
        print "$CODE_STATCHK 1>$outfn 2>$logfn" ;
        system ("$CODE_STATCHK 1>$outfn 2>$logfn");
        $statchk_exitcode = $? >> 8;
        if ( $statchk_exitcode == 0 )
        {
                if ( -s $logfn )
                {
                  $statchk_msg = "PASS";
                  print STATREPORT "Case_$TID   passed";
                }
                else
                {
                  $statchk_msg = "OK";
                  print STATREPORT "Case_$TID   succeeded";
                }
        }
        else
        {
                $statchk_msg = "FAIL";
                print STATREPORT "Case_$TID   failed";
        }
        print "[$TestDir][t$TID] Finished STATCHK, status $statchk_msg";
        chdir($my_path);
        return;                 # go read next line from configuration file
        }

# ****************************************
