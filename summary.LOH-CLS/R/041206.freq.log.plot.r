#for initial analysis of LOH-CLS data, log-y plot

rm(list=ls());

file = "INPUTFILE";

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

# output to out1
 columns = c( 2:5, 10); #columns that I will analyze.
 out1 = data.frame( matrix( nrow=n, ncol=(length(columns)+1)) ) 
 tmp.names = names(tb.m);
 names( out1 ) = tmp.names[c(1:5,10)]
 out1 = tb.m[,c(1:5,10)]

 #normalize by t=0 cell counts
 for( j in 2:length(out1[1,]) ) { #columns
  for( i in length(out1[,1]):1 ) { #row
    out1[i,j] = out1[i,j] / out1[1,j];   #normalization by t=0;
  }
 }

### plot out1 
 postscript( paste( "_", file, "ps", sep=".") );
 #cols=rainbow( length(out1[1,]) );
 cols = c("green","black","red","purple","blue");
 y.max = max( out1[, 2:length(out1[1,])], na.rm=T ) * 2;
 y.min = min( out1[, 2:length(out1[1,])], na.rm=T ) + 1e-4;
 par( yaxp=c(1e-4,20,1));
 plot(out1[,2] ~ out1$t,type='l', log='y', col=cols[1], main=file, 
   ylim=c(y.min, y.max), xlab="t",ylab="Normalized concentration");
 for ( i in 3:length( out1[1,] ) ) {
   lines( out1[,i] ~ out1$t, col=cols[i-1] );
 } 

 labels = names( out1 );
 legend( 0, 0.4 , labels[2:length(labels)], col=cols, lty=1);
 dev.off();

 q("no");

