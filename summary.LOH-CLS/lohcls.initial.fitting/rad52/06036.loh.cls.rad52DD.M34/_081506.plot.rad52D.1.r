# 
library(nlme)
rm(list=ls());

logistical.viability <- function( v, w, t ) { ret <- 1 /( 1 + ( t / v )^ w );  }

file = "060306.rad52DD1.tab";

tb= read.table( file, header=T, sep="\t", fill=T);
tb$t = tb$t / 24;
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
 row.steps = as.vector( table( tb2$t ) ); # a new trick, ha

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

# output to out,
# columns in out
 header = c("t","half.over.black","Pb","Rb","R0.5", "R0.75", "s", "g" );
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

#######calculate s , g
  for( i in length(out[,1]):1 ) { #row
    out$s[i]  = tb.m$total[i]    / tb.m$total[1]; #s
    #out[i,2] = tb.m$B0.5[i]     / tb.m$black[i];   
    out$Pb[i] = tb.m$black[i] / tb.m$total[i];
    out$g[i]  = 1 - 2 * out$Pb[i];
  }

 #postscript("060306.rad52D.1.ps")
 plot( out$s ~ out$t , type='l', main= file, col="blue");
 lines( out$Pb ~ out$t, col="black");
 labels = c("viability","black");
 ltypes = c(1,1);
 legend( (max(out$t)*0.7 ), 0.8, labels, col=c("blue", "black"), lty = ltypes);
 

######## calclulate T0.5 w_s
 t= out$t
 s= out$s
 fm.log.s <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 5, w = 4 )  );
 fm.log.s # this the half life T1/2 w_s


######## calclulate T_g and w_g
 g = out$g
 t = out$t
 g[6] = g[5] /2  + g[4] /2 ; ########## change this !!!!!!!!!
 g[11] = NA; t[11]=NA;
 g[9] =NA; t[9]=NA;

 s = s[ ! is.na(t)];
 t = t[ ! is.na(t)];
 g = g[ ! is.na(g)];

 g.min = min(g);
 g.max = max(g);

# fm.log.g <- gnls( g ~ (g.max - g.min) / ( 1 + ( t / v )^ w ) + g.min, start = list( v = 5, w = 4 )  );

# one = 1;
# fm.log.g <- gnls( g ~ (one - 0.49) / ( 1 + ( t / v )^ w ) + 0.49, params= , start = list( v = 5, w = 10 )  );
# how to use params to control gnls??? how about other functions?

# fm.log.g <- gnls( g ~ (1 - 0.49) / ( 1 + ( t / v )^ w ) + 0.49,  start = list( v = 5, w = 10 )  );

 fm.log.g # this the half life T_g w_g

 plot( s ~ t); points(g ~ t, col="red");










######i am here

quit("no");

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

