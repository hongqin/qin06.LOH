# output basal loh and max loh

rm(list=ls());
#postscript("_060306.loh.plot.rad52.1.ps",width=8, height=8);

file = "060306.rad52DD2.tab";

tb= read.table( file, header=T, sep="\t");
tb2 = tb;
labels = names( tb );

#normalize tb2
 for( j in 5:13 ) {
   for( i in  1: length(tb2[,1]) ) {
     if ( is.na(tb2[i,j] ) ) { tb2[i,j] = 0; }
     tb2[i,j] = tb2[i,j] * tb2[i,2] * tb2[i,3];
   }
 }

#generate row indice for averaging
 # row.num = c( 2, 3, 5, 6, 8, 9:15 );
 row.steps = as.vector( table( tb2$t ) );

#generate tb.m
 n.row = length( row.steps );
 col.labels = c( labels[ c(1, 5:13) ], "total" );
 #col.end = 11;

 tb.m = data.frame( matrix( nrow=n.row, ncol= length(col.labels) ) )  # mean values
 names( tb.m ) = col.labels; 

 upper.row = 0; #up pointer
 lower.row = 0; #low pointer 
 for( i in 1:n.row ) {   # for( i in 1: length(tb.m[,1]) ) {
   upper.row = lower.row + 1;
   lower.row = upper.row + row.steps[ i ] - 1 ;

   tb.m[i,1]  = tb2$t[ upper.row ]
   for( j in 2: ( length(col.labels) - 1 ) ) {
     # tb.m[i, j] = mean( tb2[ c(3*i-2, 3*i-1, 3*i), j+3] )
     tb.m[i, j] = mean( tb2[ upper.row : lower.row, j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:( length(col.labels) - 1 )], na.rm=T );
 }

######i am here
# output to out,
# columns in out
 header = c("t","half.over.black","Pb","R0.5", "R0.75", "s" );
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

 #calculate the ratios and s
  for( i in length(out[,1]):1 ) { #row
    out$s[i] = tb.m$total[i]    / tb.m$total[1]; #s
    #out[i,2] = tb.m$B0.5[i]     / tb.m$black[i];   
    out[i,3] = tb.m$black[i] / tb.m$total[i];
  }

 #postscript("060306.rad52D.1.ps")
 plot( out$s ~ out$t , type='l', main= file);
 lines( out$Pb ~ out$t, col="red");

#quit("no");

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

