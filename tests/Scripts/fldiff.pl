#The source script should be preprocessed, and one line indicating
#the location of perl (for example #!/usr/bin/perl) should be added, at its top
#
# Object: compare 2 output files from ABINIT line by line with arithmetic
# comparisons of floating point substrings
#
# Copyright (C) 1999-2020 ABINIT group (LSi,XG)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#
# Usage: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label]
#
# The first character of each line in both files indicates the mode of comparison;
# these first characters MUST be identical in corresponding lines.
# By default, a floating point comparison with a tolerance of 1.01d-10 is done
# on all floating point strings and a character comparison is done on all other
# strings, disregarding multiple spaces.
# In order for two numbers to be declared different, BOTH the absolute difference
# and the relative difference (difference divided by the sum of absolute values)
# must be bigger than the tolerance.
#
# Some special characters at the beginning of lines require a different handling:
# -	mark lines as same regardless to their content (i. e. ignore lines)
#       (can be be present in the 2 files or not, but must begin with -)
# +	mark lines as different regardless to their content
# ,	handle as + with -include option and as - with -ignore option
# P	handle as + with -includeP option and as - with -ignoreP option
# %	floating point comparisons are done with a tolerance of 1.01e-2
# ;	floating point comparisons are done irrespective of signs
# :	ignore floating point numbers and do a characters comparison
# .	do a characters comparison, but do not count this line in the Summary
#
# Both files should have the same number of non - starting lines.

# With -context option, save character strings for context and print it
# with line number when floating difference is found.
#
# The -ignore and -include options affects the treatment of the ',' 
# special character in the first column (see above) 
#
# The -ignoreP and -includeP options affects the treatment of the 'P'
# special character in the first column (see above)
#
# If -ridiculous   is specified, the default tolerance is set to 1.01e-2 
# If -easy   is specified, the default tolerance is set to 1.01e-5 
# If -medium is specified, the default tolerance is set to 1.01e-8 
# These modifications do not apply to the tolerance determined by the
# '%',and '.' first-column special signs
#
# If "label" is specified, it will appear in the last, summary line.
#
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator
$debug = 0;
# analyze options
$argptr = 0;
$ignore = 0;
$ignoreP = 0;
$include = 0;
$context = 0;
$ridiculous = 0;
$easy = 0;
$medium = 0;
while (1) {
  $ldc1 = substr($ARGV[$argptr],0,1);	# leading character
  last if ($ldc1 ne '-');	# end of options
  if ($ARGV[$argptr] eq '-ignore') {
    $ignore = 1;
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-include') {
    $include = 1;
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-ignoreP') {
    $ignoreP = 1;
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-includeP') {
    $includeP = 1;
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-context') {
    $context = 1;
    $mincontext = 4;	# minimum length for a context change
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-ridiculous') {
    $ridiculous = 1;
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-easy') {
    $easy = 1;
    $argptr ++;
    }
  elsif ($ARGV[$argptr] eq '-medium') {
    $medium = 1;
    $argptr ++;
    }
  else {
    last;	# probably file id beginning with dash
    }
  }
if ($include + $ignore >= 2) {
  print STDERR 'Error: -ignore and -include are mutually exclusive';
  print STDERR "Usage is: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label] \n";
  &summary(12);
  }
if ($includeP + $ignoreP >= 2) {
  print STDERR 'Error: -ignoreP and -includeP are mutually exclusive';
  print STDERR "Usage is: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label] \n";
  &summary(12);
  }
if ($ridiculous + $easy + $medium >= 2) {
  print STDERR 'Error: -ridiculous, -easy and -medium are mutually exclusive';
  print STDERR "Usage is: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label] \n";
  &summary(12);
  }
# set up default tolerance
$defaulttol= 1.01e-10;
$defaulttol= 1.01e-8 if ($medium == 1);        # medium value
$defaulttol= 1.01e-5 if ($easy == 1);          # easy value
$defaulttol= 1.01e-2 if ($ridiculous == 1);    # ridiculous value
$maxAbsDif = 0.0e0;
$maxRelDif = 0.0e0;
# counting all the differences independent of $tolerance
$maxAbsDif2 = 0.0e0;
$maxRelDif2 = 0.0e0;
# check file identifications
$in1 = $ARGV[$argptr];
$in2 = $ARGV[$argptr+1];
$label = $ARGV[$argptr+2];
if ($in1 eq '' || $in2 eq '') {
  print STDERR "Missing input file(s)";
  print STDERR "Usage is: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label] \n";
  &summary(12);
  }
$rc = open(FILEIN1,"<$in1");
if ($rc eq '') {
  print STDERR "Unable to open input file $in1\n";
  &summary(16);
  }
$rc = open(FILEIN2,"<$in2");
if ($rc eq '') {
  print STDERR "Unable to open input file $in2\n";
  &summary(16);
  }
# read both files a first time and check consistency
$lnm1 = 0; $ign1 = 0;
$lnm2 = 0; $ign2 = 0;
while (1) {
  while ($lin1 = <FILEIN1>) {	# read next line from file 1
    $lnm1 ++;		# bump line counter
    $ldc1 = substr($lin1,0,1);	# leading character
    $ldc1 = ' ' if ($ldc1 eq "\n" || $ldc1 eq "\r");	# handle null lines as blank
    $ldc1 = ' ' if ($ldc1 eq "_");	# treat lines beginning with "_" as normal
    last if ($ldc1 ne '-' && $ldc1 ne 'P');	# read until not to be ignored line
    $ign1 ++;		# bump ignored lines counter
    }
  while ($lin2 = <FILEIN2>) {	# read next line from file 2
    $lnm2 ++;		# bump line counter
    $ldc2 = substr($lin2,0,1);	# leading character
    $ldc2 = ' ' if ($ldc2 eq "\n" || $ldc2 eq "\r");	# handle null lines as blank
    $ldc2 = ' ' if ($ldc2 eq "_");	# treat lines beginning with "_" as normal
    last if ($ldc2 ne '-' && $ldc2 ne 'P');	# read until not to be ignored line
    $ign2 ++;		# bump ignored lines counter
    }
  if ($lin1 eq '' || $lin2 eq '') {	# end of any file ?
    $lct1 = $lnm1 - $ign1;	# count of significant lines
    $lct2 = $lnm2 - $ign2;
# end of file should occur on both files simultanously
# significant lines should correspond 1 by 1
    if ($lin1 ne $lin2 || $lct1 != $lct2) {
      print STDERR "The diff analysis cannot be done: the number of lines to be analysed differ.";
      print STDERR "File $in1: $lnm1 lines, $ign1 ignored";
      print STDERR "File $in2: $lnm2 lines, $ign2 ignored\n";
      &summary(32);
      }
    last;
    }
  if ($ldc1 ne $ldc2) {
    print STDERR "The diff analysis cannot be pursued: the leading characters differ.";
    print STDERR "File $in1, line $lnm1, $ign1 ignored, character:$ldc1";
    print STDERR "File $in2, line $lnm2, $ign2 ignored, character:$ldc2\n";
    &summary(36);
    }
  }
# close files and reread them comparing data
close(FILEIN1);
close(FILEIN2);
$rc = open(FILEIN1,"<$in1");
if ($rc eq '') {
  print STDERR "Unable to open input file $in1\n";
  &summary(16);
  }
$rc = open(FILEIN2,"<$in2");
if ($rc eq '') {
  print STDERR "Unable to open input file $in2\n";
  &summary(16);
  }
#
$lnm1 = 0;
$lnm2 = 0;
$diffcnt = 0;
$contextprv = '';
while (1) {
  $diffcntfl = 0;
  while ($lin1 = <FILEIN1>) {	# read next line from file 1
    $lnm1 ++;		# bump line counter
    $ldc1 = substr($lin1,0,1);	# leading character
    $ldc1 = ' ' if ($ldc1 eq "\n" || $ldc1 eq "\r");	# handle null lines as blank
    last if ($ldc1 ne '-' && $ldc1 ne 'P');	# read until not to be ignored line
    }
  while ($lin2 = <FILEIN2>) {	# read next line from file 2
    $lnm2 ++;		# bump line counter
    $ldc2 = substr($lin2,0,1);	# leading character
    $ldc2 = ' ' if ($ldc2 eq "\n" || $ldc2 eq "\r");	# handle null lines as blank
    last if ($ldc2 ne '-' && $ldc2 ne 'P');	# read until not to be ignored line
    }
  last if($lin1 eq '' && $lin2 eq ''); # end of both files
  if ($ldc1 ne $ldc2 && $ldc1 ne '_' && $ldc2 ne '_') {		# never be too carefull
    print STDERR "Internal error at lines $lnm1 $lnm2";
    &summary(96);
    }
# YP: ignore differences caused by MPI extra output
# FIXME: this is probably quite dirty :-s
  $lin1 =~ s/ by node([ ]+)([0-9]+)//;
  $lin2 =~ s/ by node([ ]+)([0-9]+)//;
  $lin1 =~ s/_P-[0-9]{4}_(EIG)/_$1/;
  $lin2 =~ s/_P-[0-9]{4}_(EIG)/_$1/;
# DC: in XML file, add a '.' sign at start of line for the timeInfo element.
  if ((index $lin1, "timeInfo") ge 0) {
    $lin1 = '.'.$lin1;
    $ldc1 = '.'
  }
  if ((index $lin2, "timeInfo") ge 0) {
    $lin2 = '.'.$lin2;
    $ldc2 = '.'
  }
  $same = 1;		# assume lines are similar
# + command character
  if($ldc1 eq '+') {
    $same = 0;		# handle lines as different
    }
# , command character
  elsif ($ldc1 eq ',' && $ignore+$include > 0) {
    $same = 1 if ($ignore == 1);	# handle as same and ignore if -ignore
    $same = 0 if ($include == 1);	# handle as different if -include
    }
# P command character
  elsif ($ldc1 eq 'P' && $ignoreP+$includeP > 0) {
    $same = 1 if ($ignoreP == 1);        # handle as same and ignore if -ignoreP
    $same = 0 if ($includeP == 1);       # handle as different if -includeP
    }
# _ command character
  elsif($ldc1 eq '_' || $ldc2 eq '_') {
    $same = 1;		# handle lines as same
    }
# other command: compare lines field by field, ignoring imbedded blanks
  else {
    @field1 = split(' ',substr($lin1,1));	# parse line into array
    @field2 = split(' ',substr($lin2,1));
    if ($#field1 != $#field2) {		# Error: different field counts
      print 'DBG diffldct',$lnm1,$#field1,$#field2 if($debug >= 1);
      $same = 0;
      $diffcnt ++;
      }
    else {
      if ($ldc1 eq ':') {	# no floating point tolerance - characters
        $tolerance = 0;
        }
     elsif ($ldc1 eq '.') {    # no floating point tolerance - characters, but do not count this line in the Summary
        $tolerance = 1;
        }
      elsif ($ldc1 eq '-' || $field1[0] eq '' ) { 
        $tolerance = 1;
        }
      elsif ($ldc1 eq '%') {    # floating point very high tolerance
        $tolerance = 1.01e-2;
        }
      else {			# floating point default tolerance
        $tolerance = $defaulttol;
        }
      $absoluteval=0 ;
      $absoluteval=1 if($ldc1 eq ';') ;
      print 'DBG toler',$lnm1,$tolerance if($debug >= 1);
      $contextcur = '';
      NEXTFIELD:
      for ($i = 0; $i <= $#field1; $i++) {	# compare corresponding fields
        if ($tolerance == 0) {
          if ($field1[$i] ne $field2[$i]) {	# character comparison
            print 'DBG difchar1',$lnm1,$i if($debug >= 1);
            $same = 0;
            $contextprv = '';
            $diffcnt ++;
            last;		# Error: characters fields are different
            }
          }
        elsif ($tolerance == 1) {             # special treatment : do not increment diffcnt
          if ($field1[$i] ne $field2[$i]) {     # character comparison
            $same = 0;
            $contextprv = '';
            last;               # Error: characters fields are different
            }
          }
        else {
          $_ = $field1[$i];
          $remain2 = $field2[$i];
          while ($_ ne '') {		# scan corresponding fields for floating point numbers
            &strfltrem;		# call subroutine to break field 1 into pieces
            $nafstr1 = $nafstr;
            $fltnum1 = $fltnum;
            $remain1 = $remainder;
            $_ = $remain2;
            &strfltrem;		# call subroutine to break field 2 into pieces
            if ($nafstr1 ne $nafstr) {	# compare 1st piece on both strings
# Error: characters pieces are different
              print 'DBG difchar2',$lnm1,$i if($debug >= 1);
              $same = 0;
              $diffcnt ++;
              $contextprv = '';
              last NEXTFIELD;
              }
            else {
              $contextcur = $contextcur.' '.$nafstr1 if ($nafstr1 ne '');
              }
            if ($fltnum1 ne '' && $fltnum ne '') {
# 2 corresponding floating point numbers have been found
# Might have to compute their absolute value
              if($absoluteval eq 1) {
                $fltnum1 = -$fltnum1 if ($fltnum1 < 0);
                $fltnum  = -$fltnum  if ($fltnum  < 0);
                }
# Compute difference
              $difflt =  $fltnum1-$fltnum;	
# abs function is not recognized by all versions of perl; do it inline
              $difflt = -$difflt if ($difflt < 0);
# compute sum of absolute values
              $absnum=$fltnum ;
              $absnum = -$absnum if ($absnum < 0);
              $absnum1=$fltnum1 ;
              $absnum1 = -$absnum1 if ($absnum1 < 0);
              $sumflt = $absnum1 + $absnum;
# counting all differences independent of $tolerance (GAF)
	      if($difflt > 0.0e0) {
#		  print '---BEGIN Line1=',$lnm1,'Line2=',$lnm2,'Field=',$i;
#		  print '$maxAbsDif2',$maxAbsDif2;
#		  print '$maxRelDif2',$maxRelDif2;
		  if ($difflt > $maxAbsDif2) {
		      $maxAbsDif2 = $difflt;
		  }
		  if ($sumflt > 0.0e0) {
#		      print 'DIFFLT',$difflt,'=',$fltnum1,'-',$fltnum;
#		      print 'SUMFLT',$sumflt;
#		      print '$difflt / $sumflt', $difflt,' / ',$sumflt;
		      if ($difflt*$sumflt > 0.0e0) {
			  $diffltrel2 = $difflt / $sumflt;
#			  print 'Setting: $diffltrel2',$diffltrel2
		      }
		      if ($diffltrel2 > $maxRelDif2) {
#			  print 'Line1=',$lnm1,'Line2=',$lnm2,'Field=',$i,'DIFFTREL2',$difftrel2;
#			  print 'Line1=',$lnm1,'Line2=',$lnm2,'Field=',$i,'MAXRELDIF2',$maxRelDIF2;
			  $maxRelDif2 = $diffltrel2;
		      }
		  }
#		  print '$maxAbsDif2',$maxAbsDif2;
#		  print '$maxRelDif2',$maxRelDif2;
#		  print '---END Line1=',$lnm1,'Line2=',$lnm2,'Field=',$i;
	      }
# counting only differences inside tolerance
              if ( $difflt > $tolerance && $difflt > $tolerance * $sumflt ) {
# Error: both tolerance criteria exceeded
                print 'DBG flttoler',$lnm1,$i if($debug >= 1);
                $same = 0;
#jmb                $diffcnt ++;
  		$diffcntfl ++;
		if ($difflt > $maxAbsDif) {
                  $maxAbsDif = $difflt;
		  $linAbsDif = $lnm1;
		}
                if ($sumflt > 0.0e0) {
                  $diffltrel = $difflt / $sumflt;
		  if ($diffltrel > $maxRelDif) {
                    $maxRelDif = $diffltrel;
		    $linRelDif = $lnm1;
		  }
                }
#jmb                last NEXTFIELD;
                }
              else {
                $len = length($remainder);
# drop current context if remainder is significant
                $contextcur = '' if ($len > $mincontext);
                }
              }
            elsif ($fltnum1 ne '' ^ $fltnum ne '') {
# Error: only 1 floating point number has been discovered
              print 'DBG diffloat',$lnm1,$i if($debug >= 1);
              $same = 0;
              $diffcnt ++;
              last NEXTFIELD;
              }
            $_ = $remain1;	# prepare for next piece
            $remain2 = $remainder;
            }	# end for while ($_)
          if ($remain1 ne $remain2) {
# Error: premature end of field 1
            print 'DBG endofild',$lnm1,$i if($debug >= 1);
            $same = 0;
            $diffcnt ++;
            last;
            }
          }
        }	# end for ($i
	if ($diffcntfl > 0) { $diffcnt ++; };
      }
    }
#
# print lines when different
  if ($same == 0) {
    if ($context == 1) {
      print $lnm1,$contextprv;
      }
    else {
      print $lnm1;
      }
    chop $lin1;
    print '<',$lin1;
    chop $lin2;
    print '>',$lin2;
    }
  $len = length($contextcur);
  $contextprv = $contextcur if ($context ==  1 && $len > $mincontext);
  }	# end while(1)
close(FILEIN1);
close(FILEIN2);
&summary(0);
#
sub strfltrem {
# this subroutine will try to split a string into non-floating, floating and remainder substrings
# on input: $_ = string to be parsed; contents is undefined on output
# on output: $nafstr = non-floating point substring beginning $_
#	$fltnum = floating point number if found; null string otherwise
#	$remainder = remainder substring that may be null
  local ($hit,$fixed,$expon);
  $nafstr = '';		# initialize first non-floating field
  $fltnum = '';
  $remainder = '';
  while (1) {
# split string into non-numeric and sign-dot-digit substrings
    $hit = m/([^0-9\-\+\.]*)(.*$)/;
    die "internal error 1" if (! $hit);	# assert
    $nafstr .= $1;		# concatene non-numeric substring
    print 'NAF1:',$1,$nafstr if ($debug >= 2);
    return if ($2 eq '');		# whole string is pure numeric i.e. not floating
    print "REM1:",$2 if ($debug >= 2);
    $_ = $2;
# try to split substring into fixed point and non-numeric
    $hit = m/(^[\-\+]?[0-9]*\.[0-9]+)(.*$)/;
    if (! $hit) {
# the string started with digit or sign but no fixed point number was found
# decimal point is either missing, or not followed by mantissa
      $hit = m/(.*)(^0-9\.\+\-$)/;	# split string at character that is no digit/sign/dot
      if (! $hit) {
        $nafstr .= $_;		# concatene pure numeric substring
        print 'NAF2:',$_,$nafstr if ($debug >= 2);
        return;
        }
      $nafstr .= $1;		# concatene non-numeric substring
      print 'NAF3:',$_,$nafstr if ($debug >= 2);
      $_ = $2;		# drop characters up to non-numeric
      return if ($_ eq '');		# whole string is not floating
      print "REM3:",$2 if ($debug >= 2);
      next;		# continue with remainder
      }
# a fixed point numeric string has been discovered
    print 'NAF4:',$nafstr if ($nafstr ne '' && $debug >= 2);
    $fixed = $1;
    print 'FIX:',$1 if ($debug >= 2);
    $_ = $2;
    $remainder = $2;
    print "REM4:",$2 if ($debug >= 2);
    $hit = m/(^[DdEe][\+\-]?[0-9]+)(.*)/;	# search possible exponent
    if (! $hit) {
      $fltnum = $1 +0.e0;		# convert character string to floating point number
      return;
      }
    $expon = $1;
    $expon =~ s/^[Dd]/E/;		# D not recognized by perl; change to E
    $fixed .= $expon;		# concatene mantissa and exponent
    print "FLT:",$fixed if ($debug >= 2);
    $fltnum = $fixed +0.e0;		# convert character string to floating point number
    $remainder = $2;
    print "REM5:",$2 if ($debug >= 2);
    return;
    }
  }
#
sub summary {
  local ($exitcode) = @_;
  if ($maxAbsDif2>0){
      $IncrMaxAbsDif2=$maxAbsDif2+0.001*10**int(log($maxAbsDif2)/log(10));
  }
  else
  {
      $IncrMaxAbsDif2=$maxAbsDif2
  }
  if ($maxRelDif2>0){
      $IncrMaxRelDif2=$maxRelDif2+0.001*10**int(log($maxRelDif2)/log(10));
  }
  else
  {
      $IncrMaxRelDif2=$maxRelDif2
  }

#  print STDERR "$label tolnlines= $diffcntfl tolabs= $maxAbsDif2 tolrel= $maxRelDif2";
  #printf STDERR "$label tolnlines= %2d tolabs= %.3e tolrel= %.3e\n",$diffcnt,$IncrMaxAbsDif2,$IncrMaxRelDif2;
# this routine will print a short summary of differences and exit with the code
# note that most shells only recognize integer numbers as exit codes
  if ($exitcode != 0) {
    print "Summary $label : fatal error, see above message";
    exit $exitcode;
    }
  if ($diffcnt == 0) {
    print "Summary $label : no significant difference has been found";
    exit 0;
    }
#
#  print "Summary $label : different lines= $diffcnt , max discrepancies absolute= $maxAbsDif (l.$linAbsDif), relative= $maxRelDif (l.$linRelDif)";
  printf "Summary $label : different lines= %d , max abs_diff= %.3e (l.%d), max rel_diff= %.3e (l.%d)", $diffcnt,$maxAbsDif,$linAbsDif,$maxRelDif,$linRelDif;
  exit 4;
  }
