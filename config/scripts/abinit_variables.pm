#
# PERL module to extract variables from a structured type
#

use warnings;
use strict;

sub extract_variables
{
	my ($search_file,$search_data) = @_;

	open(INP,"<$search_file")
	or die "Could not open $search_file for reading";

	my $type = "";
	my $cont = 0;
	my $wait = 1;

	my ($var,$grp,$i,@tmp1,@tmp2,%vars);

	while ( <INP> )
	{
		if ( ($wait == 0) && (/end type/) )
		{
			#print "\nEND DATASET_TYPE\n\n";
			$wait = 1;
		}

		if ( $wait == 1 )
		{
			s/^\s+//;
			if ( /^type $search_data/ )
			{
				#print "BEGIN DATASET_TYPE\n\n";
				$wait = 0;
			}
		}
		else
		{
			chomp();

			if ( $cont == 0 )
			{
				$type = "";
			}

			if ( /pointer/ )
			{
				s/,\s+pointer.*!/ ::/;
			}
			else
			{
				s/!.*//;
			}
			if ( /::/ )
			{
				@tmp1 = split(/ :: /,$_);
				$type =  $tmp1[0];
				$type =~ s/ //g;
				$_ = $tmp1[1];
			}
			if ( /\&$/ )
			{
				$cont = 1;
			}
			else
			{
				if ( $_ ne "" )
				{
					$cont = 0;
				}
			}

			s/\&//g;
			s/^\s+//;
			s/\s+$//;
			@tmp1 = split(/,/,$_);
			$grp  = 0;
			$_ = "";

			foreach $var (@tmp1)
			{
				if ( ($var =~ /\(/) && !($var =~ /\)/) )
				{
					$grp = 1;
				}

				if ( $_ eq "" )
				{
					$_ .= "$var";
				}
				else
				{
					$_ .= ",$var";
				}

				if ( ($grp == 1) && /\)/ )
				{
					$grp = 0;
				}

				if ( ($grp == 0) && ($type ne "") && ($_ ne "") )
				{
					$vars{$_} = $type;
					#print "\t<$_>: $vars{$_}\n";
					$_ = "";
				}
			}
		}
	}

	close(INP);

	return %vars;
}

1;
