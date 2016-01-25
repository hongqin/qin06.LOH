rm(list=ls());

library(nlme)
file = "050706.rad52D.2.tab";

logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) {
   (b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T ))^2 * T^w ); 
}

genome.integrity <- function(b.max, b.min, T, w, t) {
  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max);
}


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
 header = c("t","Pb","s", "g","m","Rb","R0.5", "L"  );   ### L = rate(1/2) / rate(black) 
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

#######calculate s , g
#  for( i in length(out[,1]):1 ) { #row
#    out$s[i]    = tb.m$total[i]    / tb.m$total[1]; #s
#    out$R0.5[i] = tb.m$B0.5[i] /tb.m$total[i];   
#    out$Pb[i]   = tb.m$black[i] / tb.m$total[i];
#    out$g[i]    = 1 - 2 * out$Pb[i];
#  }
 out$s      = tb.m$total / tb.m$total[1]
 out$Pb     = tb.m$black / tb.m$total;
 out$g      = 1 - 2 * out$Pb;

 out$R0.5.raw   = tb.m$B0.5 /tb.m$total; #needs to be adjusted by g.e


 ### plot of the raw data 
 plot( out$s ~ out$t , type='l', main= file, col="blue");
 lines( out$Pb ~ out$t, col="black");
 labels = c("viability","black");
 ltypes = c(1,1);
 legend( (max(out$t)*0.7 ), 0.8, labels, col=c("blue", "black"), lty = ltypes);

####### remove outliers#######################
 # out[6,] = NA;
 # out = out[ (! is.na(out$t) ), ]

######## calclulate T0.5 w_s
 t= out$t
 s= out$s
 fm.s <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 5, w = 4 )  );
 fm.s # this the half life T1/2 w_s

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 fit.m = logistical.mortality ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( out$s ~ out$t );
 lines( fit.s ~ t, col="blue");
 lines( fit.m ~ t, col="brown");

######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)
 b = out$Pb
 t = out$t;
 b.min = min( b ); b.min  # 0.09297052
 b.max = max( b ); b.max  # 0.2163743

 #####   Pr[B>=b] = b.max - (b.max - b.min) /( 1 + ( t / T )^ w ) ########formula for cumulative Pr(b)
 fm.b = gnls( b ~ 0.2163743 - (0.2163743 - 0.09297052) / (1 + (t/T.g)^w ), start=list( T.g=2, w=2) );
 fm.b
 plot( fm.b );

 t = seq(0, max(out$t),by=0.1);
 fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 
 plot( out$Pb ~ out$t );
 lines( fit.b ~ t, col="blue");

 ###################overlay s, m, Pb, an R0.5, half/full, mortality
 #postscript("090806.060306rad52D.1.ps", width=8, height=8)
 par(mar=c(5,4,4,4)+0.1);
 
 plot( out$s ~ out$t, col="blue", xlab="t(days)",ylab="percentage" );
 title ("rad52DD 050706 #1");                                        ############change here
 lines( fit.s ~ t, col="blue",lty=2);
 points( out$Pb ~ out$t );
 lines( fit.b ~ t, lty=2);

######calculate Rb = d(Pb)/dt
 # derivative.black <- function(b.max, b.min, T, w, t) {  (b.max - b.min) * w * t ^(w -1) / (1 + (t / T ))^2; }
 out$dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );
 out$g.e = genome.integrity(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );

 out$Rb     = out$dPb / out$g.e;                        ### rate of becoming blacks
 out$R0.5   = out$R0.5.raw / out$g.e;                   #### 090706 change

 out$L = out$R0.5 / out$Rb

 ## lines ( out$Rb ~ out$t, col="brown", lty=2) 
 lines ( out$R0.5 ~ out$t, col="green", lty=2); 
 points ( out$R0.5 ~ out$t, col="green");
 
 par(new=T);
 plot( out$L ~ out$t, col="red",type='l', axes=F, xlab="", ylab="" ); 
 axis(4, at=pretty(range(out$L[2:length(out$L)]) ) )
 mtext( "ratio", 4, 2);

 #add labels here
 cols  =c("blue","black","green","red");
 labels=c("viability","P(b)","R(1/2)","half/full)")
 ltypes=c(2,2,2,1) 
 legend( 5 , 8 ,labels, col=cols, lty=ltypes);

 dev.off();

# quit("yes");


