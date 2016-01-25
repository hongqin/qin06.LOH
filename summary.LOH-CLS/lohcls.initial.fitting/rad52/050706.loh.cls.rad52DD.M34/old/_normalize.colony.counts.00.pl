#!/usr/bin/perl
# kw parse_raw_lifespan_data.00.pl 
# v0.0 060206Fri

# Input are csv files from Excel, which contains empty space and various number of time points.
# I will average the rows from the same time points and replace empty space with 0

print "not implemented. bye."; exit;

use strict; use warnings; use Getopt::Long;

if (! $ARGV[0]) { _help(); exit(1); }

my $in_file  = undef;	
my $out_file = undef;
my $NA_output_flag = 0; 	
my $debug = 0;

GetOptions ('i|in_file=s' =>  \$in_file,    # required, no path in the file name
            'o|out_file=s' => \$out_file,   # optional
	    'na|NA_output=i' => \$NA_output_flag,  	# default 0
            'd|debug=i' => \$debug);

if ( ! $in_file )  { _help(); exit(1); }
if ( ! $out_file )  { 
 # $out_file = "$in_file.tab";   
 my @names = split( /\//, $in_file );
 my @els = split( /\./, $names[$#names] );
 pop @els; pop @els; 
 # print "@els";
 $out_file = join( ".", @els );
 $out_file .= ".new";  ### 062404 changes
}

#variables
my (@lines, $line ) = ();
my %time_points = ();

open (IN, "<$in_file" ) ;
while ( $line = <IN> ) {
  if (! ( $line =~ /^\s*$/ ) ) {
    push ( @lines, $line );  
  }
}
close (IN);

my $header = shift @lines; 

if ($debug) { print "\@lines=[@lines]\n"; }
  
#set up %time_points
for ( my $i=0; $i<= $#lines; $i++) {
 my @els = split (/\t/, $lines[$i]);
 %time_points{ $els[0] } = $i . "\t";
}

my @sorted_lifespans = _sort (@lifespans);

close (OUT);
#
# DONE
#

sub _help {
print "\n  $0 \n
GetOptions ('i|in_file=s' =>  \$in_file,    # required, no path in the file name
            'o|out_file=s' => \$out_file,   # optional
	    'na|NA_output=i' => \$NA_output_flag,  	# default 0
            'd|debug=i' => \$debug);
"
}

sub _sort {
 my @num = ( @_ );
 for ( my $i=0; $i<=($#num-1); $i++ ) {
   for ( my $j= $i + 1; $j<=$#num; $j++ ) {
     if ( $num[$i] > $num[$j] ) {
         my $tmp = $num[$j]; $num[$j] = $num[$i]; $num[$i]=$tmp;
     }
   }
 }
 return @num;
}
