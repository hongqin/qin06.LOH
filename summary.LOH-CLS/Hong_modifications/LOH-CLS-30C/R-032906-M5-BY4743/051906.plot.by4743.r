# output basal loh and max loh

rm(list=ls());
#postscript("05176.M34.ps",width=8, height=8);

file = "BY4743.032906.tab";

tb= read.table( file, header=T, sep="\t");
row.num = c(1:30, 33,33,33, 36,36,36, 39,39,39, 40:42, 45,45,45, 47,47,48, 50,50,51, 52:78)	
tb2 = tb[row.num, ];
labels = names( tb );

#normalize tb2
 for( j in 5:length(tb2[1,]) ) {
   for( i in  1: length(tb2[,1]) ) {
     if ( is.na( tb2[i,j] ) ) { 
	tb2[i,j] = 0;
     } else {
        tb2[i,j] = tb2[i,j] * tb2[i,2] * tb2[i,3];
     }
   }
 }

#generate tb.m
 n=length(tb2$t) /3 
 col.end = 10;

 tb.m = data.frame( matrix( nrow=n, ncol= col.end ) )  # mean values
 names( tb.m ) = c( labels[ c(1, 5:12) ], "total" );

 for( i in 1: length(tb.m[,1]) ) {
   tb.m[i,1]  = tb2$t[3*i -1]
   for( j in 2: 9 ) {
     tb.m[i, j] = mean( tb2[ c(3*i-2, 3*i-1, 3*i), j+3] )
     #tb.sd[i,j] = sd( tb2[ c(2*i-1, 2*i), j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:9], na.rm=T );
 }

# output to out,
# columns in out
 header = c("t","half.over.black","Pb","R0.5", "R0.75", "s" );
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

 #calculate the ratios and s
  for( i in length(out[,1]):1 ) { #row
    out$s[i] = tb.m$total[i]    / tb.m$total[1]; #s
    out[i,2] = tb.m$B0.5[i]     / tb.m$black[i];   
  }


 #calculate the risks (in contrast to percentages)
  for( i in length(out[,1]):1 ) { #row
    out$Pb[i]    = tb.m$black[i]  / tb.m$total[i];  # full blacks
    out$R0.5[i]  = tb.m$B0.5[i]   / ( tb.m$total[i] - 2 * tb.m$black[i] );  # half blacks
  }

###generate plots
 tokens = unlist( strsplit(file, "/") );
 shortfile = tokens[1];
 #postscript( paste( "_", shortfile, "loh.042506.ps", sep=".") );

 par(mar=c(4,4,4,4));
 #plot the ratio
 y.max = max( out[,2], na.rm=T ) * 1.2;
 y.min = min( out[,2], na.rm=T ) ;
 cols = c("blue","black","green");
 if ( y.max == Inf ) { y.max = 7; };
 plot(out[,2] ~ out$t,type='l', col=cols[1],
   main= paste( shortfile, " ", "Half vs Full Blacks" ),
   ylim=c(y.min, y.max), xlab="t (hours)",ylab="Half/Full Ratio");

  #overlay viability plot. 
  par(new=T);
  plot( out$s ~ out$t,type='l',xlab="",ylab="",axes=F,lty=2); 	

 #plot the percentage and risks
 y.max = max( out[,3:4], na.rm=T ) * 1.2;
 y.min = min( out[,3:4], na.rm=T ) ;
 par(new=T);
 plot(out$Pb ~ out$t, type='l',  col=cols[2],
   ylim=c(y.min, y.max), xlab="",ylab="",axes=F);
 lines( out[,4] ~ out$t, col=cols[3] );
 axis(4, at=pretty(c(y.min, y.max)) );
 mtext("Percentage", side=4, line=3);

 #legend
 labels = c("Half/Full", "P(full)", "R(1/2)", "viability");
 ltypes = c(1,1,1,2);
 legend( 0, y.max * 0.8, labels, col=c(cols,"black"), lty = ltypes);

 #dev.off();
#} ###file loop####

#dev.off();
# q("no");

