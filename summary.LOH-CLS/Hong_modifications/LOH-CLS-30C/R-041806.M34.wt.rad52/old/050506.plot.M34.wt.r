# output basal loh and max loh

rm(list=ls());
#postscript("05056.loh.plot.m34.wt.ps",width=8, height=8);

#file = "rad52.041806.tab";
file = "M34.041806.tab";

tb= read.table( file, header=T, sep="\t");
row.num = c(1:12,rep(15,3), rep(18,3), rep(21,3), 25,26,26, 28:30, rep(33,3) );  ###
row.num = 1 : length( tb[,1]);
tb2 = tb[row.num, ];
labels = names( tb );

#adjust by dilutions
for( i in 1:length(tb2[,1]) ) {
 for (j in 5:13) {
   tb2[i,j] = tb2[i,j] * tb2[i,2] * tb2[i,3];
 }
}


#generate tb.m
 n=length(tb2$t) / 3
 col.end = 10;

 tb.m = data.frame( matrix( nrow=n, ncol= col.end ) )  # mean values
# names( tb.m ) = c( labels[ c(1, 5:12) ], "total" );
 names( tb.m ) = labels[ c(1, 5:13) ]

 for( i in 1:n ) {
   tb.m[i,1]  = tb2$t[3*i -2]
   #tb.sd[i,1] = tb.m[i, 1]
   for( j in 2:col.end ) {
     tb.m[i, j] = mean( tb2[ c(3*i-2, 3*i-1, 3*i), j+3] )
     #tb.sd[i,j] = sd( tb2[ c(3*i-2, 3*i-1, 3*i), j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:col.end], na.rm=T );
 }

# output to out,
# columns in out:
 header = c("t","half.over.black","Pb","R0.5", "R0.75", "s" );
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out[,1] = tb.m[,1]; # "t"


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
 shortfile = tokens[2];
 #postscript( paste( "_", shortfile, "loh.042506.ps", sep=".") );

 par(mar=c(4,4,4,4));
 #plot the ratio
 y.max = max( out425[,2], na.rm=T ) * 1.2;
 y.min = min( out425[,2], na.rm=T ) ;
 cols = c("blue","black","green");
 if ( y.max == Inf ) { y.max = 7; };
 plot(out425[,2] ~ out425$t,type='l', col=cols[1],
   main= paste( shortfile, " ", "Half vs Full Blacks" ),
   ylim=c(y.min, y.max), xlab="t (hours)",ylab="Half/Full Ratio");

  #overlay viability plot. 
  par(new=T);
  plot( out425$s ~ out425$t,type='l',xlab="",ylab="",axes=F,lty=2); 	

 #plot the percentage and risks
 y.max = max( out425[,3:4], na.rm=T ) * 1.2;
 y.min = min( out425[,3:4], na.rm=T ) ;
 par(new=T);
 plot(out425$Pb ~ out425$t, type='l',  col=cols[2],
   ylim=c(y.min, y.max), xlab="",ylab="",axes=F);
 lines( out425[,4] ~ out425$t, col=cols[3] );
 axis(4, at=pretty(c(y.min, y.max)) );
 mtext("Percentage", side=4, line=3);

 #legend
 labels = c("Half/Full", "P(full)", "R(1/2)", "viability");
 ltypes = c(1,1,1,2);
 legend( 0, y.max * 0.8, labels, col=c(cols,"black"), lty = ltypes);

 #dev.off();
} ###file loop####

dev.off();
# q("no");

