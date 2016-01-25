# remove outlier?

rm(list=ls());

library(nlme)
file = "060306.rad52DD1.tab";

logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

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
 #old# header = c("t","half.over.black","Pb","Rb","R0.5", "R0.75", "s", "g" );
 header = c("t","Pb","s", "g","Rb","R0.5", "R0.75",  );
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

####### remove outliers#######################
 out[6,] = NA;
 out = out[ (! is.na(out$t) ), ]

######## calclulate T0.5 w_s
 t= out$t
 s= out$s
 fm.s <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 5, w = 4 )  );
 fm.s # this the half life T1/2 w_s

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( out$s ~ out$t );
 lines( fit.s ~ t, col="blue");

######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)
 b = out$Pb
 t = out$t;
 b.min = min( b ); b.min  #0.1688
 b.max = max( b ); b.max  #0.25504

 fm.b = gnls( b ~ 0.25504 - (0.25504 - 0.1688) / (1 + (t/T.g)^w ), start=list( T.g=2, w=5) );
 fm.b
 plot( fm.b );

 t = seq(0, max(out$t),by=0.1);
 fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 
 plot( out$Pb ~ out$t );
 lines( fit.b ~ t, col="blue");



# quit("yes");


