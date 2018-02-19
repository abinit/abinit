# The purpose of this script is to get a summary
# of the calculated ab initio energies (and components)
# forces, stresses ... in a whole set of files
# Output will be written in the file "report"
#
# Copyright (C) 1993-2018 ABINIT group (XG,LSi)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .

# Written by X. Gonze, January 8, 1993, this script has been translated to Perl
# by L. Sindic.
#
# Last modifications: 2005/10/12   Y. Pouillon
#
# Usage :
# unix C-shell: make-cpu-report file ... >& log_file 
# unix bash:	make-cpu-report file ... >log_file 2>&1
# Windows DOS box: [perl] make-cpu-report.pl file ... >log_file
#
# Of course, you can use metacharacters (e.g. make-cpu-report *), in order
# to scan a very large set of files in a directory

# To set up the information wanted in the report, mention all key words in the
# keyword array ( @keyword = ... )
# Then for each key word, define the range of lines to be printed in the
# linemin and linemax arrays. These are the numbers of the first and last line
# relatively to the line where the key word stands, being counted as 1.
# Note that only the first occurrence of a keyword will give an output in the
# report, not the subsequent ones.
# For example, if you want the atomic positions, when two atoms are present,
# from an 'output' file, you need to get two lines, the first being identified
# by the keyword ' xred ', you would write:
#	@keyword = (' xred ');
#	@linemin = (  1   );
#	@linemax = (  2   );
#
# Here is a more elaborate example with multiple keywords :
#	@keyword=('getcut',' xred ','dE','Ewald');
#	@linemin=(   1    ,   1    , 1  ,  1   );
#	@linemax=(   2    ,   2    , 3  ,  1   );
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator

# Here set up the information; 
@keyword = ('getcut','Arith','Band','Detailed','Overall' );
@linemin = (   1    ,   1   ,   1  ,    1     ,  1 );
@linemax = (   2    ,   1   ,   1  ,   12     ,  1 );

# Here begins the code
# ***********************************************************
$FNsummary = 'report';		# default file name for summary
# Count the number of files (shell-expanded arguments list)
$numfiles = $#ARGV + 1;
print "The following files will be analysed :";
for ($ifile = 0; $ifile < $numfiles; $ifile++) {
	print " $ARGV[$ifile]";
	}
print " Number of files = $numfiles\n";

# Issue the list of keywords and their range
$numberkw = $#keyword + 1;	# count of keywords
print "The following keywords will be searched, with range:";
for ($indexkw = 0; $indexkw < $numberkw; $indexkw++) {
	if ($linemin[$indexkw] < 1 || $linemax[$indexkw] < $linemin[$indexkw]) {
		print "Invalid range $linemin[$indexkw] $linemax[$indexkw] for keyword $keyword[$indexkw]";
		exit 16;
		}
	write;		# print "keyword linemin linemax" according to format STDOUT
	}
# Clean the report file if it exists; otherwise, create it 
open(REPORT,">$FNsummary") || die "Unable to open file $FNsummary in write mode";

# Set up a loop on all the files
for ($ifile = 0; $ifile < $numfiles; $ifile++) {
# Examine each file
	$file = $ARGV[$ifile];
	print " File $file is now treated ";
	open(FILEI,"<$file") || die "Unable to open input file $file";
	$startprt = 0; $stopprt = 0;	# initialize flag-counters
# Read one line until end of file
	while (<FILEI>) {
		chop $_;	# drop NL
# Search line for each keyword	
		for ($indexkw = 0; $indexkw < $numberkw; $indexkw++) {
			$ipos = index($_,$keyword[$indexkw]);	# scan line
			next if ($ipos < 0);	# keyword not found; try next one
			if ($startprt <= 0) {	# not printing for now ?
				$startprt = $linemin[$indexkw];	# first line to print
				$stopprt = $linemax[$indexkw];	# last line to print
				}
			else {
# another keyword has been recognized while processing the lines in a range of a
# previous hit; merge both ranges :
				$startprt = $linemin[$indexkw] if ($startprt > $linemin[$indexkw]);
				$stopprt = $linemax[$indexkw] if ($stopprt < $linemax[$indexkw]);
				}
			}	# end of loop on the keywords
		if ($startprt == 1) {	# is current line within range ?
			print REPORT $_ ;	# add the information to the energies file
			if ($startprt == $stopprt) {	# last line in range ?
				$startprt = 0;		# reset flag-counters
				$stopprt = 0;	
				}
			}
# if the current line is still short of 1st one in range, bump down its relative number:
		$startprt-- if ($startprt > 1);
		$stopprt-- if ($stopprt > 1);	# bump down relative number of last line
		}		# end of reading loop
# Add one separation line in report between the different files:
	print REPORT " ";	# add a blank line
	close(FILEI);
	}			# end of loop on the files
close(REPORT);
exit;
# Format to print keyword and ranges:
format STDOUT =
@<<<<<<<<<<<<<<<<<<< @>>> @>>>
$keyword[$indexkw],$linemin[$indexkw],$linemax[$indexkw] 
.
