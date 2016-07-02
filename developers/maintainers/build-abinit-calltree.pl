#!/usr/bin/perl

use strict;
use warnings;

sub get_rank
{
	($_) = @_;

	if ( /^[0-9]/ )
	{
		s/^([0-9]*).*/$1/;
		s/^0([0-9])/$1/;
	}
	else
	{
		$_ = -1;
	}

	return $_;
}

if ( $#ARGV == -1 )
{
	print "Usage: build-abinit-calltree directory\n";
	exit 0;
}

my $root_dir = $ARGV[0];

my @source_dirs;
my %points_from;
my %points_to;

opendir(SRC,"$root_dir/src");
while ( $_ = readdir(SRC) )
{
	my $current_dir = $_;

	if ( /^[0-9]/ )
	{
		print "Reading $current_dir (rank ".get_rank($current_dir).")...";

		push(@source_dirs,$current_dir);

		opendir(USE,"$root_dir/src/$current_dir");
		while ( $_ = readdir(USE) )
		{
			my $current_file = $_;

			if ( /\.F90$/ )
			{
				open(INP,"<$root_dir/src/$current_dir/$current_file");
				while ( <INP> )
				{
					if ( (/use interfaces_/) && !(/use interfaces_$current_dir/) )
					{
						chomp();
						s/.*interfaces_//;
						s/,.*//;
						if ( $points_to{$current_dir} )
						{
							if ( !(join(' ',@{$points_to{$current_dir}}) =~ $_) )
							{
								push(@{$points_to{$current_dir}},$_);
							}
						}
						else
						{
							$points_to{$current_dir} = [$_];
						}
						if ( $points_from{$_} )
						{
							if ( !(join(' ',@{$points_to{$current_dir}}) =~ $_) )
							{
								push(@{$points_from{$_}},$current_dir);
							}
						}
						else
						{
							$points_from{$_} = [$current_dir];
						}
					}
				}
				close(INP);
			}
		}
		closedir(USE);

		print "done.\n";
	}
}
closedir(SRC);

my ($src,$dep);

open(OUT,">abinit-calltree.dot");
foreach $src (sort @source_dirs)
{
	print OUT "digraph ABINIT_$src {\n   ratio = fill;\n"
                . "   rankdir = LR;   size = \"6.3,10.2\";\n"
	        . "   label = \"Dependencies on and of $src\";\n";

	if ( $points_from{$src} )
	{
		foreach $dep (sort @{$points_from{$src}})
		{
			if ( get_rank($dep) <= get_rank($src) )
			{
				print OUT "   src_$dep [color=red,fontcolor=red,style=bold];\n";
			} 
			print OUT "   src_$dep -> src_$src;\n";
		}
	}

	if ( $points_to{$src} )
	{
		foreach $dep (sort @{$points_to{$src}})
		{
			if ( get_rank($dep) >= get_rank($src) )
			{
				print OUT "   src_$dep [color=red,fontcolor=red,style=bold];\n";
			} 
			print OUT "   src_$src -> src_$dep;\n";
		}
	}

	print OUT "}\n\n";
}
close(OUT);

system("dot -Tps -o abinit-calltree.ps abinit-calltree.dot");
system("ps2pdf abinit-calltree.ps");
system("pdfnup --nup 2x1 abinit-calltree.pdf");

#vim:ts=4
