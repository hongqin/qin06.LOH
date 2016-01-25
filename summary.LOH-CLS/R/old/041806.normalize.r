#for initial analysis of LOH-CLS data, regular plot.

rm(list=ls());

file = "INPUTFILE";
#file = "M28.022006.txt";

tb = read.table( file, header=T, sep="\t");

tb2 = tb;
labels = names( tb );

# normalize the data, to tb2
 for( j in 5: 12 ) {tb2[,j] = tb2[,j] * tb2[,2] * tb2[,3]  }

# reorganize the data by days, to tb.m
 n = length(tb2$t) / 3
 col.end=9;

 tb.m = data.frame( matrix( nrow=n, ncol=9) )  # mean values
 names( tb.m ) = labels[ c(1, 5:12) ];
 #tb.sd = tb.m   # standard deviations

 for( i in 1:n ) {
   tb.m[i,1]  = tb2$t[3*i -2]
   #tb.sd[i,1] = tb.m[i, 1]
   for( j in 2: 9 ) {
     tb.m[i, j] = mean( tb2[ c(3*i-2, 3*i-1, 3*i), j+3] )
     #tb.sd[i,j] = sd( tb2[ c(3*i-2, 3*i-1, 3*i), j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:9], na.rm=T );
 }

# todo: remove NA and zeros from tb.m

#output tb.m
 outfile = paste( "normalized.data/norm.", file, sep="");
 write.table( tb.m, outfile, row.names=F, sep="\t", quote=F);

# q("no");

