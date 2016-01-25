#! /usr/bin/perl 
# overly R0.5/P(black), P(black), R0.5, and viability

print "aborted, switch to R.\n";
exit 0;

BEGIN { unshift(@INC,"/home/hong_man/lib/perl/", "/home/hqin/lib/perl");   }

use strict; use warnings; use Util;

my $file = "_normalized.data.file";
my $script = "R/042506.loh.plot.r";

my (@files, @strains, $line) = ();
open(IN, "<$file"); 
while( $line=<IN> ) {
 if ($line !~ /^\s+$/) {
  chomp $line;
  my @tokens = split (/\t/, $line);
  push( @files,   $tokens[0]);
  push( @strains, $tokens[1]);
 }	
}
	close(IN); 
print @files; print "\n";
print @strains;print "\n";

#generate working R script
open(IN, "<$script"); my @lines=<IN>; close(IN);
my $super_script = join( "\n", @lines);

for(my $i=0; $i<=$#files; $i++) {
   my $current_file = "normalized.data/$files[$i]";
   my $strain = $strains[ $i ];
   my $new_super_script = $super_script;
   $new_super_script =~ s/INPUTFILE/$current_file/g;
   my $timestamp = get_short_time_stamp_US();
   my $tmp_script = "/tmp/_$strain.lohrisk.$timestamp.r";
   open(OUT,">$tmp_script"); print OUT $new_super_script; close (OUT);

   #my @tokens = split( /\./, $fl);
   my $cmnd="R --no-save < $tmp_script";
   print "\$current_file=[$current_file]\n$cmnd\n";
   system($cmnd);
   sleep 1;
}

exit 0;
	
