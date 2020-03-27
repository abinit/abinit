
# This script runs one abinit built-in test in background mode while
# periodically displaying the status file till the end of the process.

# Copyright (C) 1999-2020 ABINIT group (LSi,XG)
# This file is distributed under the terms of the
# GNU General Public License, see ~ABINIT/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~ABINIT/Infos/contributors .
#
# Usage :
# This script is intended to be called from Makefile (unix) or make.bat
# (DOS/Windows). Execute "make help" for details. It can also be called
# directly:
# unix shell: run-basic-tests machine_name [ 1 | 2 | 3 | 4 | 5 | 6] srcDir
# Windows DOS box: perl run-basic-tests.pl machine_name [ 1 | 2 | 3 | 4 | 5 | 6] srcDir
# MacOS (likely...) : perl -f run-basic-tests.pl machine_name [ 1 | 2 | 3 | 4 | 5 | 6] srcDir
#
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator
#
$DELAY = 3;		# delay in seconds between 2 verifications of the STATUS file
$NCHECKS = 10;		# max number of times to check
$MAXLOOP = '50';	# maximum count for STATUS file checking loop
%PspFiles=('in_fast','01h.pspgth',
	 'in_v1','70yb.pspnc',
	 'in_v5','01h.pspgth:04be.pspgth',
	 'in_bigdft','01h.pspgth',
	 'in_etsf_io','20ca.paw',
	 'in_libxc','83bi.psphgh',
	 'in_wannier90','31ga.pspnc:33as.pspnc');
#
$debug = 0;		# Debug level ???
if ($ARGV[1] eq '') {
  print 'Missing argument';
	print "Usage is: run-basic-tests machine_name [ in_fast | in_v1 | in_v5 | in_bigdft | in_etsf_io | in_libxc | in_wannier90 ] srcDir";
	exit 16;
	}
$MACH = $ARGV[0];
$TestN = $ARGV[1];
$srcDir = $ARGV[2];

$TMPfile = 'Test_in.tmp';	# multi usage temporary file
$UNIXSLASH = '/';	# unix subdirectory delimitor in file paths (octal 057)
$DOSSLASH = '\\';	# subdirectory delimitor in DOS file paths (\ escaped)
$MACSLASH = ':';  # subdirectory delimitor in MAC file paths
$MULTI = '*';
# fetch some environment variables and try unix 'uname' command (Bourne shell)
$OSname = $ENV{'OS'};		# OS name on Windows systems
$OStype = $ENV{'OSTYPE'};	# defined on some systems as OPENVMS, mingw/msys
# Use internal variable of Perl $^O
if ($^O eq 'MacOS') {
  $OStype = 'MacOS'; $rc = 0;
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
	$OStype = 'DOS';	# set OStype
# since unlink <file*> fails on some versions of perl for DOS, let's use del
	$SLASH = $DOSSLASH;	# subdirectory delimitor for file paths
	$COPY_COMMAND = 'copy';         # DOS copy file command
	$ERASE_COMMAND = 'del /q';	# DOS delete file command
	$TYPE_COMMAND = 'type';		# DOS type file command
# Define the DOS conventions to call the dotest perl script by means of the
# "system" function
# NOTE: perl.exe is a DOS module that MUST be accessible through the DOS PATH
	$SUFXstyle = 'DOS';	# use DOS-style suffixes for binaries
	$PLSUFX = '.pl';	# DOS suffix for perl script
	$PRLPFX = 'start /min perl ';	# DOS command for background
	$BGSUFX = '';		# no special DOS suffix for new task
	$ENDtest = $MAXLOOP;	# loop limit in case dotest fails early
	print "File test$TestN.end will be checked to test end of abinit execution";
	print "The wait loop of this script is limited to $ENDtest cycles";
	}
    elsif ( $UNkernel eq 'CYGWIN_NT' || $UNkernel eq 'MINGW32_NT' || $UNkernel eq 'MINGW64_NT' ) {
# unix-like environment under Windows
	$OStype = $UNkernel; 	# set OStype for cygwin/mingw
# subdirectory delimitor for file paths is DOS \ under cygwin and standard /
# under other unices
	$SLASH = $UNkernel eq 'CYGWIN_NT' ? $DOSSLASH : $UNIXSLASH;
	$COPY_COMMAND = 'cp -p';  # unix copy file command
	$ERASE_COMMAND = 'rm -f';	# unix delete file command
	$TYPE_COMMAND = 'cat';		# unix type file command
	$SUFXstyle = 'unix';	# use unix-style suffixes for binaries,.
	$PLSUFX = '.pl';	# no special suffix for perl script under unix
	$PRLPFX = 'perl ';	# perl path ??? should use $(PERL)
	$BGSUFX = '&';		# special suffix for background task
	$ENDtest = 'ps';	# use unix ps command to check end of test
	}
    else {
	print "unrecognized Windows Subsystem =$OStype=$UNkernel=";
	exit (99);
	}
    }
# if not Windows NT check other environment variables and uname output
elsif ($OStype eq 'OPENVMS') {
    $SLASH = $UNIXSLASH;		# subdirectory delimitor in file paths
    $COPY_COMMAND = 'copy';     # OpenVMS copy file command
    $ERASE_COMMAND = 'delete ';	# OpenVMS delete file command
    $TYPE_COMMAND = 'type ';		# OpenVMS type file command
# for perl under normal unix systems:
    $SUFXstyle = 'vms';	# use unix-style suffixes for binaries, ...
    $PLSUFX = '.pl';		# no special suffix for perl script under unix
    $PRLPFX = 'perl ';		# perl path defined in first line of script
    $BGSUFX = '';		# special suffix for background task
    $ENDtest = 'vms';	# use gnv ps command to check end of test
    }
# MacOS section
elsif ($OStype eq 'MacOS') {
    $SLASH = $MACSLASH;             # subdirectory delimitor in file paths
    $COPY_COMMAND = 'Duplicate -c'; # copy file command
    $ERASE_COMMAND = 'Delete -y';   # delete file command
    $SUFXstyle = $OStype;        # use unix commands but MAC-style suffixes
    $PLSUFX = '.pl'; 	       # suffix for perl script
    }
# normal unix/linux section
elsif ($OSname eq '' && $UNkernel ne '') {
    $OStype = $UNkernel;	# set OStype for *ux
    $SLASH = $UNIXSLASH;	# subdirectory delimitor in file paths
    $COPY_COMMAND = '\cp ';      # copy file command
    $ERASE_COMMAND = '\rm -f';	# unix delete file command
    $TYPE_COMMAND = 'cat';	# unix type file command
# for perl under normal unix systems:
    $SUFXstyle = 'unix';	# use unix-style suffixes for binaries, ...
    $PLSUFX = '.pl';		# XG060708 : .pl is needed ... ?? no special suffix for perl script under unix
    $PRLPFX = 'perl ' ;		# XG060708 : perl is needed ... ?? perl path defined in first line of script
    $BGSUFX = '&';		# special suffix for background task
    $ENDtest = 'ps';	# use unix ps command to check end of test
    }
else {
    print "unrecognized Operating System -$OSname=$OStype=$UNkernel-";
    exit (99);
    }
# make sure "files" file exist
$inputfile = &transpath("$srcDir/built-in/Input/test$TestN.in");
print "$OStype $SLASH $SUFXstyle $inputfile" if ($debug >=2);
if (! -e &transpath("$inputfile") ) {
	print "Invalid test number: $TestN";
	print "Expected input file $inputfile not found";
	print "Usage is: run-basic-tests machine_name [ in_fast | in_v1 | in_v5 | in_bigdft | in_etsf_io | in_libxc | in_wannier90 ] srcDir ";
	exit (8);
	}
#
print "Built-in test $TestN will be run through run-basic-dotest$PLSUFX script" if ($debug >= 1);
# set a date flag to make a new directory today #
($sec,$min,$hour,$mday,$ymon,$yyear,$wday,$yday,$isdst)=localtime(time);
$ymon++;        # ymon was 0-11
$yyear +=1900;  # yyear was relative to 1900
$YYYYMMDD = sprintf("%4.4d",$yyear).sprintf("%2.2d",$ymon).sprintf("%2.2d",$mday);
$WORK_DIR = 'built-in';
if (! -e $WORK_DIR || ! -d $WORK_DIR) {
  mkdir ($WORK_DIR,0755);         # Mode 0755 ignored under DOS-Windows
  }
else {
  print "Do not create directory, $WORK_DIR already exists" if ($debug >= 1);
  }
print ('cd built-in') if ($debug >= 1);
chdir ('built-in');
$WORK_DIR = 'tmp-'.$MACH.'_'.$OStype.'_'.$YYYYMMDD;
if (! -e $WORK_DIR || ! -d $WORK_DIR) {
  mkdir ($WORK_DIR,0755);         # Mode 0755 ignored under DOS-Windows
  }
else {
  print "Do not create directory, $WORK_DIR already exists" if ($debug >= 1);
  }
print "cd $WORK_DIR" if ($debug >= 1);
chdir ("$WORK_DIR");

# suppress old files
unlink("test$TestN.end") if (-e "test$TestN.end");
unlink("test$TestN.log") if (-e "test$TestN.log");
unlink("test$TestN.out") if (-e "test$TestN.out");
unlink("test$TestN".'.o_WFK') if (-e "test$TestN".'.o_WFK');
if ($OStype eq 'OPENVMS') {
  $statfile = "test$TestN".'_STATUS.dat';	# name of STATUS file
  }
else {
  $statfile = "test$TestN".'_STATUS';	# name of STATUS file
  }
unlink($statfile) if (-e $statfile);
# Write input files for test
system("$COPY_COMMAND ".&transpath("$srcDir/built-in/Input/test$TestN.in")." test$TestN.in");
if ($TestN eq 'in_wannier90') {
  system("$COPY_COMMAND ".&transpath("$srcDir/built-in/Input/testin_wannier90o_w90.win")." testin_wannier90o_w90.win");
  }
# Write abinit.files for test
$rc = open (FILES,">test$TestN.files");
print FILES "test$TestN.in";
print FILES "test$TestN.out";
print FILES "test$TestN".'i';
print FILES "test$TestN".'o';
print FILES "test$TestN";
@PsP = split(':', $PspFiles{"$TestN"} );
print "Psp files: @PsP" if ($debug >=1);

if( $UNkernel eq 'MINGW64_NT') {
   print FILES &transpath("/cygwin$srcDir/Psps_for_tests/$PsP[0]");
   print FILES &transpath("/cygwin$srcDir/Psps_for_tests/$PsP[1]") if ($PsP[1] ne '');
} else {
   print FILES &transpath("$srcDir/Psps_for_tests/$PsP[0]");
   print FILES &transpath("$srcDir/Psps_for_tests/$PsP[1]") if ($PsP[1] ne '');
}

print FILES &transpath("$srcDir/Psps_for_tests/$PsP[0]");
print FILES &transpath("$srcDir/Psps_for_tests/$PsP[1]") if ($PsP[1] ne '');
close (FILES);
# abinit will be started through a dotest perl script
# this is necessary under DOS/Windows because of file redirections
$cmd = $PRLPFX.&transpath("$srcDir/Scripts/run-basic-dotest$PLSUFX")." $TestN $ENDtest $SLASH $SUFXstyle $BGSUFX";
print $cmd if ($debug >= 2);
$rc = system($cmd);
die "Error $rc starting run-basic-dotest" if ($rc != 0);
# loop waiting end of abinit execution
$ichk = $NCHECKS;
while ($ichk > 0) {
	sleep($DELAY);		# take a few seconds rest

	if ($ENDtest eq 'ps') {
# check abinit process using unix ps
		$cmd = 'ps | grep abinit | grep -v grep > '.$TMPfile;
		print $cmd if ($debug >= 2);
		$rc = system($cmd);
		last if (-z $TMPfile);	# no abinit process was found
		}
	elsif ($ENDtest eq 'vms') {
# always end for VMS
    last;
		}
	else {
# check presence of file "test$Test0N.end"  (like test1.end)
		last if (-e "test$TestN.end");	# dotest has completed
		$ENDtest --;	# decrement loop counter
		if ($ENDtest <= 0) {
			print 'Wait loop limit reached, run-basic-dotest not yet finished';
			exit;		# loop limit reached
			}
		}
	$cmd = "$TYPE_COMMAND $statfile";
	system($cmd) if (-e $statfile);	# type status file if it exists
	$ichk--;
	}
#
$cmd = "$TYPE_COMMAND $statfile";
system($cmd) if (-e $statfile);	# type status file if it exists
unlink ("$TMPfile") if (-e $TMPfile);
exit (0);

# ****************************************
sub transpath {
  local($path) = @_;
#
# purpose: translate unix-like path of file to DOS-like according to $SLASH
# argument:
#       $path = path to be translated
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
