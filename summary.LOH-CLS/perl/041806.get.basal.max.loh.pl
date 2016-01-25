#! /usr/bin/perl 

BEGIN { unshift(@INC,"/home/hong_man/lib/perl/", "/home/hqin/lib/perl");   }

use strict; use warnings; use Util;

my $filedir = "loh.risk";
my $debug = 1;

my (@files, @strains, $line) = ();
my $time_stamp = get_short_time_stamp_US();

opendir DIR, $filedir;
@files = readdir DIR;
closedir DIR;

my @loh_files = grep (/^loh.norm/, @files);
if ($debug) { 
 foreach my $fl (@loh_files) { print "$fl ";} print "\n"; 
 print "There are ".($#loh_files+1)." files.\n";
}

chdir $filedir;
open (OUT, ">../_loh.risk.summary.tab");
print OUT "exprt\tPb\tR0.5\tR0.25\tR2\n";
foreach my $fl (@loh_files) {
 open (IN, "<$fl"); my @lines = <IN>;  close (IN);
 my @tokens = split( /\./, $fl ); shift @tokens; shift @tokens; pop @tokens;
 my $experiment = join( '.', @tokens);
 my $header = shift @lines;
 my ($t,$a,$b,$full, $half, $quarter, $two, @rest) = split( /[\t\n]+/, $lines[0]);
 print OUT "$experiment\t$full\t$half\t$quarter\t$two\n";
}
close(OUT);


exit 0;
	
