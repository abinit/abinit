# Check that all main subroutines follow ABINIT requirements

use strict;
use warnings;

# Find out source subdirs
my @dirs = ["10_defs"];
opendir(DIR,"./src");
while ( $_ = readdir(DIR) )
{
	if ( /^[0-^9]/ )
	{
		$dirs[++$#dirs] = "./src/$_";
	}
}
closedir(DIR);
$dirs[++$#dirs] = "98_main";

# Setup counters
my $nsrc = 0;
my $ncod = 0;
my $nchr = 0;
my $nref = 912; # Number of character* statements in version 4.3.3
my $nsub = 0;
my $nfun = 0;
my $ngto = 0;
my $ndos = 0;
my $nshl = 0;
my $nsok = 0;
my $nend = 0;
my ($dir,$src);

# Explore each subdir
foreach $dir (@dirs)
{
	# Find out source files
	my @srcs = [];
	opendir(DIR,"$dir");
	while ( $_ = readdir(DIR) )
	{
		if ( /\.F90$/ )
		{
			$srcs[++$#srcs] = $_;
		}
	}
	closedir(DIR);

	# Look meticulously inside each source file
	foreach $src (@srcs)
	{
		# Setup internal variables
		my $name = "";
		my $type = "";
		my %args;
		my @lst;
		my $arg;

		my $lvl = 0;

		# Read file
		open(SRC,"<$dir/$src");
		#print "Entering $dir/$src\n";
		while ( <SRC> )
		{
			# Count code lines
			if ( ! (/^!/ || /^$/) )
			{
				$ncod += 1;
			}

			# Count character* statements
			if ( /character\*/ || /character\(\*\)/ )
			{
				$nchr += 1;
			}

			# Clean-up input line
			chomp();
			s/^\s+//;
			s/\s+$//;
			s/\'.*//;
			s/\".*//;
			s/\!.*//;
			s/\#.*//;

			# Parsing level 2 : inside subprograms (must be first!)
			if ( $lvl == 2 )
			{
				# Look for intent statements
				if ( /::/ )
				{
					foreach $arg (keys %args)
					{
						if ( /intent/ && /$arg/ )
						{
							$args{$arg} = 1;
						}
					}
				}

				# Look for goto statements
				if ( /^goto/i )
				{
					$ngto += 1;
					print "$dir/$src($name): contains goto statements\n";
					print STDERR "WARNING: goto statement found inside $name in $dir/$src\n";
				}

				# Look for shared loops
				if ( /^do\s/i || /do$/i )
				{
					$ndos += 1;
				}
				if ( /^do\s+[0-9]+/i )
				{
					$nshl += 1;
				}

				# Catch subprogram end
				if ( /^end$/i || /^end (subroutine|function)/i || /^interface/i || /^contains/i )
				{
					# finish with 'end' is not a good programming practice
					if ( /^end$/i )
					{
						print STDERR "NOTICE: $type $name in $dir/$src ends with 'end'\n";
						$nend += 1;
					}

					# Subprograms should be named at their end
					if ( /^end (subroutine|function)$/i )
					{
						$nend += 1;
					}

					# Count all arguments and intented ones
					my $ttl = 0;
					my $int = 0;

					foreach $arg (keys %args)
					{
						$ttl += 1;
						$int += $args{$arg};
						delete $args{$arg};
					}

					if ( $int == $ttl )
					{
						$nsok += 1;
					}
					else
					{
						print "$dir/$src($name): $int/$ttl argument(s) with intent\n";
					}
					#print "\tExiting $name\n";

					$lvl = 0;
				}
			}

			# Detect subprogram starts
			if ( (/^subroutine/i || /function/i) && ! m/^end/i )
			{
				# Parsing level should be 0 : outside subprograms
				if ( $lvl != 0 )
				{
					print STDERR "WARNING: Bad parsing level $lvl after $name in $dir/$src\n";
				}

				# Determine and count subprogram type
				if ( /subroutine/i )
				{
					$type = "subroutine";
					$nsub += 1;
				}
				if ( /function/i )
				{
					$type = "function";
					$nfun += 1;
				}

				# If subprogram has arguments, go to parsing level 1
				if ( m/\(/i )
				{
					$lvl = 1;
				}

				# Get the name of the subprogram
				$name = $_;
				$name =~ s/.*(subroutine|function)//i;
				$name =~ s/\(.*//;
				$name =~ s/\s+//;

				#print "\tEntering $name\n";
			}

			#  Parsing level 1 : identifying arguments
			if ( $lvl == 1 )
			{
				# Closing parenthesis means that we are now
				# inside the subprogram
				if ( /\)/ )
				{
					$lvl = 2;
				}

				# Clean-up input line
				s/.*\(//;
				s/\).*//;
				s/\&//g;
				s/\s+//;
				s/\!.*//;

				# Get arguments and hash them
				@lst = split(/,/,$_);
				#print "$src($name): args = <@lst>\n";

				foreach $arg (@lst)
				{
					$args{$arg} = 0;
				}
			}
		}

		close(SRC);
		$nsrc += 1;
	}
}

# Final report
print "-- \n";
printf "%3d source files explored (%d lines of code).\n",$nsrc,$ncod;
printf "%3d subprograms identified (%d subroutines, %d functions).\n",
	$nsub+$nfun,$nsub,$nfun;
if ( $ngto != 0 )
{
	printf "%3d goto statements found.  *** SHOULD BE 0 !!! ***\n",$ngto;
}
printf "%3d 'character*' statements found (%.1f%% done).\n",
	$nchr,(1.0-$nchr*1.0/$nref)*100.0;
printf "%3d shared loops over %d have been found (%.1f%% done).\n",
	$nshl,$ndos,(1.0-$nshl*1.0/($ndos))*100;
printf "%3d subprograms end without name (%.1f%% done).\n",
	$nend,(1.0-$nend*1.0/($nsub+$nfun))*100;
printf "%3d subprograms do not have intented parameters (%.1f%% done).\n\n",
	$nsub+$nfun-$nsok,($nsok*1.0/($nsub+$nfun))*100;

