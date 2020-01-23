#
# Object: translate Fortran 77 source code to free format Fortran 90
# Usage: FT77to90 [-dp] [-lf] F77SOURCE
#	-dp will replace the double precision statements using kind
# 	-lf will replace labelled format statements by character strings
#
# Copyright (C) 1998-2020 ABINIT group (LSi)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#
#
# F77SOURCE should be a Fortran 77 source file
# Translated output file name:
# if source file is 77.f-suffixed, output file suffix will be changed to 90.f
# if source file is .f77-suffixed, output file suffix will be changed to .f90.f
# if source file is .f-suffixed, output file suffix will be changed to .f90.f
# in any other case, suffix .f90.f will be appended to source file name
#
# Translation to free format actually features the following:
# 1) C, c or * in column 1 is replaced by ! for comments;
# lines starting with # in column 1, especially cpp statements, are left unchanged
# 2) line numbering in columns 73-80 is suppressed when possible
#  RESTRICTIONS:
#   a) sequence numbers that start on a comment are not suppressed
#   b) sequence numbers on format statements are not suppressed when labelled
# formats replacement is enabled (-lf option)
#  NB: Fortran 90 allows line length to be extended from 80 to 132 characters
# 3) old style relational operators (even if found in comments but not on cpp stmt)
#    .EQ. .GE. .GT. .LE. .LT. .NE.   are replaced as follows:
#     ==   >=   >    <=   <    /=  
#  RESTRICTION: relational operators that span two lines won't be changed
# 4) records continuation with non blank in column 6 is replaced by ampersands
# (esperluette in french) at end and beginning of segments;
#  RESTRICTION: cpp directives will be handled as comments, leading to possible
# compile-time errors if cpp directives stay between continued and continuation
# records !
# 5) DOUBLE PRECISION declarations are replaced using KIND (-dp option required)
# 6) labelled formats will be changed to character strings that will be
# initialized throuh parameter statements (-lf option required)
# the source file will be pre-scanned to find and save format statements,
# and to find a suitable place to insert their character strings declarations;
# During the normal scan, labelled formats will be suppressed, formats will be
# declared as parameters character strings and labelled format references
# in read and write statements will be changed accordingly
#  RESTRICTIONS to labelled format substitution:
# 	a) format not followed by ( on the same line is ignored
# 	b) print, accept and type statements are not processed
# 	c) double quotes " in a format are not processed correctly: escaping
# through replication is not performed
# 	d) a declaration statement like "implicit","integer","double",
# "parameter" or "dimension" must be found in the source file so that
# character strings format declarations will be inserted before any
# executable statement
#	e) this procedure assumes there is only one subroutine, function or
# program statement in the source file; if this is not the case, format
# declarations may be inserted in a wrong program unit
# 	f) tabs are not expanded; some lines may be misinterpreted as
# continuation of a format and suppressed
#
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator
$ix = 0;
$dblprec = 0;
$suplabfmt = 0;
# check parameters
while(substr($ARGV[$ix],0,1) eq '-') {
  if ($ARGV[$ix] eq '-dp') {
    $dblprec = 1;
    $ix ++;
    next;
    }
  elsif ($ARGV[$ix] eq '-lf') {
    $suplabfmt = 1;
    $ix ++;
    next;
    }
  else {
    print "invalid option $ARGV[$ix]";
    exit 24;
    }
  }
$in = $ARGV[$ix];
if (-r $in eq '') {
  print "Unreadable source file: $in";
  exit 12;
  }
$out = $in;
if (substr($in,-4,4) eq '77.f') {
  substr($out,-4,4) = '90.f';	# change 77.f to 90.f
  }
elsif (substr($in,-4,4) eq '.f77') {
  substr($out,-2,2) = '90.f';	# change .f77 to .f90.f
  }
elsif (substr($in,-2,2) eq '.f') {
  $out .= '90.f';	# change .f to .f90.f
  }
else {
  $out .= '.f90.f';	# append suffix
  }
if (-e $out) {
  print "Output file $out already exits";
  exit 8;
  }
$rc = open(FILEIN,"<$in");
if ($rc eq '') {
  print "Unable to open input file; error $rc";
  exit 16;
  }
#
if ($suplabfmt == 1) {
# source file will be read a first time to search labelled
# format definitions and save them in associative arrays
  $linnum = 0;          # source line number
  $fmtread = 0;         # not yet reading/saving a format
  $lblfmt = 0;          # initialize counter
  $lastdcnum = 0;       # no declaration found yet
  $lastdclin = '';
  $saveline = '';       # buffer containing previous line
  while (1) {
    if ($saveline eq '') {
      $line = <FILEIN>;
      last if ($line eq '');      # end-of-file
      $linnum ++;         # update line number
      }
    else {
      $line = $saveline;
      $saveline = '';
      }
    $col1 = substr($line,0,1);
# skip comments, cpp statements: c, C or * in column 1
    next if ( $col1 eq 'c' || $col1 eq 'C' || $col1 eq '*' || $col1 eq '#');
    $isn = substr($line,6);
# handle continuations
    if (length($line) > 6 && substr($line,5,1) ne ' ') {
      next if ($fmtread == 0);  # skip continuation if not a format
      chop $isn;
      $fmtnum ++;                       # initialize format lines count
      %fmtlines = (%fmtlines,
        join($;,$label,$fmtnum),$isn);  # save line in 2-dim array
      $fmtlinct{$label} = $fmtnum;      # update lines count
      next;
      }
    else {
      $fmtread = 0;             # reset reading/saving flag
      }
    &checkformat($line,$label,$ip);
    if ($label eq '') {         # not a labelled format
# check for declaration statement
      $col1to6 = substr($line,0,6);
      ($wd1,$wd2,$wd3,$wd4) = split(' ',$isn);
      if ( $col1to6 eq '      ' && ( $wd1 eq 'implicit' || $wd1 eq 'integer' || $wd1 eq 'double' || $wd1 eq 'parameter'
|| $wd1 eq 'dimension' || substr($wd1,0,9) eq 'character') ) {
        &nextstmt($saveline);  # skip to next statement
        $lastdcnum = $linnum;   # save line number of statement following this declaration
        $lastdclin = $saveline; # save this line
#DBG print $lastdcnum,$lastdclin;
        }
      next;
      }
    $isn = substr($line,$ip);
    if ($fmtlinct{$label} ne '') {
      print "Warning ! format",$label,"on line",$linnum,"ignored: duplicate label";
      next;
      }
    chop $isn;
    $fmtnum = 1;                # initialize format lines count
    %fmtlines = (%fmtlines,
      join($;,$label,$fmtnum),$isn);    # save line in 2-dim array
    %fmtlinct = (%fmtlinct,     # save lines count in associative array
      $label,$fmtnum);
    %fmtrefct = (%fmtrefct,     # initialize reference count
      $label,0);
    $lblfmt ++;                 # bump counter
    $fmtread = 1;               # remember reading/saving a format
    }
#
# close the source file and read it again
  close(FILEIN);
# check if format declaration statements may be inserted somewhere:
  if ($lblfmt > 0 && $lastdcnum == 0) {
    print "Error: unable to insert $lblfmt format(s); no declaration statement found";
    exit 20;
    }
  }
#
$rc = open(FILEIN,"<$in");
if ($rc eq '') {
  print "Unable to open input file; error $rc";
  exit 16;
  }
#
$rc = open(FILEOUT,">$out");
if ($rc eq '') {
  print "Unable to open output file; error $rc";
  exit 16;
  }
$linnum = 0;            # source line number
$declfmt = 0;           # format declation statements not yet inserted
$line = '';
$seqnum = '';		# initialize numbering off
$f90comnt = '!';	# FORTRAN 90 comment prefix
$cppdir = '#';		# cpp directive prefix
$prevcol1 = '';		# initialize column1 of previous record
$writepos = 0;		# write pointer at Begin of File
%relatops = (
	'\.eq\.' , '==' ,
	'\.ge\.' , '>=' ,
	'\.gt\.' , '>' ,
	'\.le\.' , '<=' ,
	'\.lt\.' , '<' ,
	'\.ne\.' , '/=' );		# define translation table for relational operators
# dots are escaped by \ for regular expression substitutions
#
while ( $_ = <FILEIN>) {	# read next line
  $saveline = $line;	# save previous line if any
  $line = $_;		# move input
  $len = length($_);
  $linnum ++;		# bump line counter
  $col1 = substr($line,0,1);
#
# 6a) labelled formats processing: insert format declaration/initialization
# statements into code after last declaration statement detected
  if ($suplabfmt == 1 && $lblfmt > 0 && $linnum == $lastdcnum) {
    if ($lastdclin ne $line) {
      print "Error: last declaration statement at line $linnum is not";
      print $lastdclin;
      exit 64;
      }
    foreach (sort keys %fmtlinct) {
      $dclfmt = "      character(*) format".$_."\n";
      $len2 = length($dclfmt);
      &putline($dclfmt,$len2,' ');	# write line to output file
      for ($i = 1;$i <= $fmtlinct{$_}; $i++) {
        $prfx = $i == 1 ? '      parameter (format'.$_.' ="' : '     &';
        $sufx = $i == $fmtlinct{$_} ? '")' : '&';
        $ix = index($fmtlines{$_,$i},'"');
        print 'Warning ! double quote found in format',$_ if ($ix >= 0);
        $dclfmt = $prfx.$fmtlines{$_,$i}.$sufx."\n";	# build subsequent lines
        $len2 = length($dclfmt);
        &putline($dclfmt,$len2,' ');	# write line to output file
        }
      }
    $declfmt = 1;               # format declarations inserted
    }
#
# 6b) labelled formats processing: skip labelled formats
  if ($suplabfmt == 1) {
    &checkformat($line,$label,$ip);
    if ($label ne '') {
#DBG  print "skipping format", $label,$linnum;
      &nextstmt($saveline);      # skip format up to next statement
      last if ($saveline eq '');	# leave loop on end-of-file
      $line = $saveline;
      $len = length($line);
      }
    }
#
# 1) handle comments
  if ($line =~ /^c|^C|^\*/) {		# comments: c, C or * in column 1 ?
    substr($line,0,1) = $f90comnt;	# substitute !
    }
  $col1 = substr($line,0,1);
#
# 2) check fixed format and suppress sequence numbers when possible
  if ($len == 81) {
    $_ = substr($line,72,8);	# extract columns 73-80
# check for an 8-digits number
    if (/\d{8}/) {
      if ( $seqnum ne '') {
        if ($seqpfx eq 'none' && $_ > $seqnum) {
#	suppress numbering if same format and increasing sequence:
          $line = substr($line,0,72)."\n";
          $len -= 8;		# reduce line length by 8
          $incr = $_ - $seqnum; # sequence increment
          $seqnum = $_;		# update sequence number
          }
        else {			# broken sequence
          $seqnum = '';		# turn numbering off
#	a new suppress sequence may start for the same line by following code:
          }
        }
      if ($seqnum eq '' && $col1 ne $f90comnt && $col1 ne $cppdir) {
# start suppressing numbers if numbering was off and line is neither a comment
# nor a cpp statement: 
        $line = substr($line,0,72)."\n";
        $len -= 8;		# reduce line length by 8
        $seqnum = $_;		# set numbering on
        $seqpfx = 'none';
        }
      }
    if (/\D{3}\d{5}/) {
# check for a 3-chars prefix followed by a 5-digits number
      $seqnum2 = substr($_,3,5);	# extract prefix and number
      $seqpfx2 = substr($_,0,3);
      if ($seqnum ne '') {
        if ( $seqpfx2 eq $seqpfx && $seqnum2 > $seqnum) {
#	suppress numbering if same format and increasing sequence:
          $line = substr($line,0,72)."\n";
          $len -= 8;		# reduce line length by 8
          $incr = $seqnum2 - $seqnum;	# sequence increment 
          $seqnum = $seqnum2;	# update sequence number
          }
        else {			# broken sequence
#	two consecutive lines with the same numbering were found once !
          $seqnum = '';		# turn numbering off
#	a new suppress sequence may start for the same line by following code:
          }
        }
      if ($seqnum eq '' && $col1 ne $f90comnt && $col1 ne $cppdir) {
# start suppressing numbers if numbering was off and line is neither a comment
# nor a cpp statement: 
        $line = substr($line,0,72)."\n";
        $len -= 8;		# reduce line length by 8
        $seqnum = $seqnum2;	# set numbering on
        $seqpfx = $seqpfx2;
        }
      }
    }
  else { $seqnum = ''; }	# non fixed format: turn numbering off
#
# 3) handle old relational operators, cpp statements excepted
  if ($col1 ne $cppdir) {
    foreach $relop (keys(%relatops)) {
      $line =~ s/$relop/$relatops{$relop}/ig;	# ignore case, global change
      }
    }
  $len = length($line);	# update length
#
# 4) handle double precision statements: they will be declared REAL with the
# KIND parameter. Since the values for this parameter are not defined by
# Fortran and may vary from one processor to another, the special intrinsic
# inquiry function KIND will be used on the value 0.0d0
# (cfr FORTRAN 95 by Martin Counihan page 40)
  if ($dblprec == 1 && $col1 ne $f90comnt && $col1 ne $cppdir) {
    $line =~ s/^\s*double\s+precision\s+/      real(kind=kind(0.0d0)) /i;
    $len = length($line);	# update length
    }
#
# 6c)labelled formats processing: check for read/write statements;
# when format is referenced by its label, replace it by character variable name
if ($suplabfmt == 1) {
  $isn = substr($line,6);
  $pntr = 6;                    # remember column pointer
# check for read/write (UNIT,FORMAT_LABEL) [iolist]
  if ( $col1 ne $f90comnt) { 	# skip comment
# look up for: "word ( remainder "
    ($left1,$remain1) = split(/\(/,$isn,2);
    ($wd1,$wd2) = split(' ',$left1);
    if (($wd1 eq 'read' || $wd1 eq 'write') && $wd2 eq '') {    # found read / write (
      $pntr += index($isn,'(') + 1;     # point to remainder
# split control information list into specifiers:
      ($spec1,$spec2,$spec3,$spec4,$spec5,$spec6,$spec7) = split(',',$remain1);
      $ix = index($spec2,')'); 
      $spec2 = substr($spec2,0,$ix) if ($ix > 0);       # drop possible )
      $len2 = length($spec2);
      $spec2 =~ tr/ //d;           # strip off blanks
      $labnum = $spec2 + 0;        # convert possible label to numeric
      $spec2 =~ tr/0123456789//d;  # does digits suppression yield null string ?
      $label = sprintf("%5.5d",$labnum);        # change to 5 digits number
      if (index($spec1,'=') < 0 && $spec2 eq '') {
        if ($fmtlinct{$label} eq '') {
          print "Error: unable to replace labelled format for $wd1 on line $linnum: undefined label $labnum";
          }
        else {
          $pntr += index($remain1,',') +1;      # point to spec2
# replace label reference by character string name:
          substr($line,$pntr,$len2) = 'format'.$label;
          $fmtrefct{$label} += 1;               # bump reference counter
          $len = length($line);		# update length
          }
        }
# check for read/write ( [specifier,...] fmt=FORMAT_LABEL, [specifier,...]) [iolist]
      else {
        $ix = index($remain1,'fmt=');
        if ($ix >= 0) {
          $pntr += $ix + 4;     # point past =
          $remain2 = substr($line,$pntr);
# truncate specifier to next , or )
          $ixvirg = index($remain2,',');
          $remain2 = substr($remain2,0,$ixvirg) if ($ixvirg > 0);
          $ixpard = index($remain2,')');
          $remain2 = substr($remain2,0,$ixpard) if ($ixpard > 0);
          $len2 = length($remain2);
          $remain2 =~ tr/ //d;           # strip off blanks
          $labnum = $remain2 + 0;        # convert possible label to numeric
          $remain2 =~ tr/0123456789//d;  # does digits suppression yield null string ?
          $label = sprintf("%5.5d",$labnum);        # change to 5 digits number
          if ($remain2 eq '') {
            if ($fmtlinct{$label} eq '') {
              print "unable to replace labelled format for $wd1 on line $linnum: undefined label $labnum";
              }
            else {
# replace label reference by character string name:
              substr($line,$pntr,$len2) = 'format'.$label;
              $len = length($line);         # update length
              $fmtrefct{$label} += 1;               # bump reference counter
              }
            }
          }
        }
      }
    }
  }
#
# 5) handle continuations
# possible mishandling of a line containing TAB character:
  if ($col1 ne $f90comnt && index($line,"\t") >= 0) {
    print "Warning ! source line $linnum containning TAB may be mishandled; line=";
    print $line;
    }
  $col6 = substr($line,5,1);	# possible continuation character
# 	WARNING ! cpp directives will be handled as a comments
# Fortran allows comments between continued and continuation lines:
# so, if previous line was not a comment but current line is, insert a blank
# at the end of previous line; this blank could be overlaid by an & if the
# current line has a continuation following the comment(s)
  if ($prevcol1 ne '' && ($prevcol1 ne $f90comnt && $prevcol1 ne $cppdir) && ($col1 eq $f90comnt || $col1 eq $cppdir)) {
    $rc = seek (FILEOUT,-1,2);  # backspace one character
    if ($rc == 0) {
      print "Error backspacing file $out";
      exit 100;
      }
    $writepos --;		# backspace
    $line = " \n".$line;	# prepend line with blank-nl
    $len += 2;			# bump length
    }
  if ($col1 ne $f90comnt && $col1 ne $cppdir && $len >6 && $col6 ne ' ') {	# process continuation
    substr($line,5,1) = '&';	# substitute & in column 6 of continuation
    if ($endnoncom ne '')  {	# funny program starting by a continuation ?
# comments lay between continued line and continuation: seek to last blank
# character of continued statement and susbtitute &
      if ($prevcol1 eq $f90comnt || $prevcol1 eq $cppdir) {
        $rc = seek (FILEOUT,$endnoncom,0);  # seek to end of last non comment
        if ($rc == 0) {
          print "Error seeking file $out to position $endnoncom";
          exit 100;
          }
        $rc = syswrite (FILEOUT,'&',1); # overlay blank with &
        if ($rc == 0) {
          print "Error $rc writing to $out file";
          exit 100;
          }
        $rc = seek (FILEOUT,0,2);  # seek to end of file
        if ($rc == 0) {
          print "Error seeking $out to end-of-file";
          exit 100;
          }
        }
      else {
# continuation follows immediately continued line: backspace one character,
# prepend current line with &-newline to continue previous line
        $rc = seek (FILEOUT,-1,2);  # backspace one character
        if ($rc == 0) {
          print "Error backspacing file $out";
          exit 100;
          }
        $writepos --;		# backspace
        $line = "&\n".$line;	# prepend line with &-nl
        $len += 2;		# bump length
        }
      }
    }
#
  &putline($line,$len,$col1);	# write buffered line
  }
close(FILEIN);
close(FILOUT);
#
# 6d) labelled formats processing: make sure work has been done and check for
# unreferenced format
if ($suplabfmt == 1) {
  foreach (keys %fmtrefct) {
    print "Warning ! labelled format $_ is unreferenced" if ($fmtrefct{$_} == 0);
    }
# make sure format declaration have not been lost:
  if ($declfmt == 0 && $lblfmt > 0) {
    print "Error: $lblfmt formats were not inserted";
    exit 20;
    }
  }
exit;
#
# check for labelled format statement
sub checkformat {
  $cfline = $_[0];
  $_[1] = '';           # no numeric label for now
  $colabl = substr($cfline,0,5);        # numeric label area
  return if ($colabl eq '     ');
# NB: '1 2' in label zone defines label 12 !
  $colabl =~ tr/ //d;           # strip off blanks
# NB: '001' or '  1' in label zone define a single label 1 !
  $labnum = $colabl + 0;        # convert possible label to numeric
  $colabl =~ tr/0123456789//d;  # does digits suppression yield null string ?
  return if ($colabl ne '');    # ignore if not
  $cfisn = substr($cfline,6);   # instruction area
  $ixformat = index($cfisn,'format');   # search keyword
  return if ($ixformat < 0);            # ignore if not a format
  if ($ixformat > 0) {
    $_ = substr($cfisn,0,$ixformat);
    return if (/\S/);           # not a format if not blank
    }
  $label = sprintf("%5.5d",$labnum);    # change to 5 digits number
  $ixparen = index($cfisn,'(',$ixformat+6);
  if ($ixparen < 0) {
    print "format",$label,"on line",$linnum,"ignored: no parenthesis";
    return;
    }
  $_[1] = $label;
  $_[2] = 6 + $ixparen;
#DBG print "format",$labnum,$label,"has been found",$cfisn;
  return;
  }
#
# process continuation lines up to next statement
sub nextstmt {
  while ( $_ = <FILEIN>) {      # read next line
    $linnum ++;                 # update line number
    $column1 = substr($_,0,1);     # comment column
    $column6 = substr($_,5,1);
# check for non comment continuation
    next if (length($_) >= 6 && $column1 ne 'c' && $column1 ne 'C' && $column1 ne '*' && $column6 ne ' ');
    $_[0] = $_;
    return;
    }
  $_[0] = '';           # return end of file
  return;
  }
#
# putline (line, length, column1)	# write a line to output file
sub putline {
  $rc = syswrite (FILEOUT,$_[0],$_[1]) ; # write buffered line
  if ($rc != $_[1]) {
    print "Error $rc writing to $out file";
    exit 100;
    }
  $prevcol1 = $_[2];    # remember comment or not
  $writepos += $_[1];    # bump write position
# save write position of last instruction that might be continued:
  $endnoncom = $writepos-1 if ($_[2] ne $f90comnt && $_[2] ne $cppdir);
  return;
  }
