
rm(list=ls());

library(nlme)
file = "030105.M34.lohcls.csv";

logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) {
   (b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w ); 
}

genome.integrity <- function(b.max, b.min, T, w, t) {
  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max);
}

t = seq(1,15, by=0.1)
T = 5; w = 10;
b.max = 1
b.min =0;
Pb = logistical.black(b.max, b.min, T, w, t ) 
plot(Pb~t)
dPb = derivative.black(b.max, b.min, T, w, t ) 
plot(dPb~t)






tb= read.table( file, header=T, sep="\t", fill=T);
tb$t = tb$t / 24;
tb2 = tb;
labels = names( tb );

tb2$B0.5[22]=1;  ###########a upper limit for B0.5

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
 for( i in 1:n.row ) {
   upper.row = lower.row + 1;
   lower.row = upper.row + row.steps[ i ] - 1 ;

   tb.m[i,1]  = tb2$t[ upper.row ]
   for( j in 2: ( length(col.labels) - 1 ) ) {
     tb.m[i, j] = mean( tb2[ upper.row : lower.row, j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:( length(col.labels) - 1 )], na.rm=T );
 }

 ### get the standard errors of whites and blacks #################### 091306 change
 upper.row = 0; #up pointer
 lower.row = 0; #low pointer 
 sd.w    =  numeric( n.row );
 sd.b    = numeric( n.row );
 sd.b05 = numeric( n.row);
 for( i in 1:n.row ) { 
   upper.row = lower.row + 1;
   lower.row = upper.row + row.steps[ i ] - 1 ;
   if ( ( lower.row - upper.row ) > 0 ) {
     sd.w[i]    = sd( tb2$white[ upper.row : lower.row] ) ;
     sd.b[i]    = sd( tb2$black[ upper.row : lower.row] ) ;
     sd.b05[i]  = sd( tb2$B0.5[ upper.row : lower.row] ) ;
   } 
 }
 tb.m$sd.w = sd.w
 tb.m$sd.b = sd.b
 tb.m$sd.b05 = sd.b05


# output to out,
# columns in out
 #old# header = c("t","half.over.black","Pb","Rb","R0.5", "R0.75", "s", "g" );
 header = c("t","Pb","s", "g","m","Rb","R0.5", "L"  );   ### L = rate(1/2) / rate(black) 
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

#######calculate s , g
 out$s      = tb.m$total / tb.m$total[1]
 out$Pb     = tb.m$black / tb.m$total;
 out$g      = 1 - 2 * out$Pb;
 out$e.s    =  out$s * (tb.m$sd.w / tb.m$white);
 out$e.b    =  out$Pb * (tb.m$sd.b / tb.m$black);

 out$R0.5.raw   = tb.m$B0.5 /tb.m$total; #needs to be adjusted by g.e
 out$e.b05  =  out$R0.5.raw * ( tb.m$sd.b05 / tb.m$B0.5)

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
 fm.s <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 5, w = 11 )  );
 fm.s # this the half life T1/2 w_s

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 fit.m = logistical.mortality ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( out$s ~ out$t );
 lines( fit.s ~ t, col="blue");
 lines( fit.m ~ t, col="brown");

 #estimate errors for s
 error.s = out$s -  logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], out$t );
 error.s = abs(error.s)
 out$e.s = ifelse( out$e.s==0, error.s, out$e.s );
 out$e.s= ifelse( out$e.s<0.01, 0.01, out$e.s );
 
######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)
 b = out$Pb
 t = out$t;
 b.min = min( b ); b.min  # 0.001262626
 b.max = max( b ); b.max  # 0.2394366

 #####   Pr[B>=b] = b.max - (b.max - b.min) /( 1 + ( t / T )^ w ) ########formula for cumulative Pr(b)
 fm.b = gnls( b ~ 0.2394366 - (0.2394366 - 0.001262626) / (1 + (t/T.g)^w ), start=list( T.g=6, w=9) );
 fm.b

 #estimate error for black
 error.b = out$Pb - logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], out$t);
 error.b = abs( error.b );
 out$e.b = ifelse( out$e.b==0, error.b, out$e.b );
 out$e.b = ifelse( out$e.b<0.01, 0.01, out$e.b);

 plot( fm.b );

 t = seq(0, max(out$t),by=0.1);
 fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 
 plot( out$Pb ~ out$t );
 lines( fit.b ~ t, col="blue");

 ###################overlay s, m, Pb, an R0.5, half/full, mortality
 postscript("091706.030105M34.ps", width=8, height=8)
 par(mar=c(5,4,4,4)+0.1);
 
 plot( out$s ~ out$t, col="blue", xlab="t(days)",ylab="percentage", ylim=c(-0.05, 1.1), pch=16 );
 title ("M34 030106");                                        ############change here
 arrows( out$t, (out$s - out$e.s), out$t, (out$s + out$e.s), length=0.1, angle=90,code=3, lty=1, col="blue" );

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 lines( fit.s ~ t, col="blue",lty=2);
 
 points( out$Pb ~ out$t, pch=15 );
 arrows( out$t, (out$Pb - out$e.b), out$t, (out$Pb + out$e.b), length=0.1, angle=90,code=3, lty=1 );
 lines( fit.b ~ t, lty=2);

 cols  =c("blue","black" );
 labels=c("viability","P(b)")
 ltypes=c(2,2,) 
 pch   =c(16,15)
 legend( 9 , 1 ,labels, col=cols, lty=ltypes, pch=pch);

 dev.off();

##################################3
 #postscript("091706.030105M34.overlay.ps", width=8, height=8)
 par(mar=c(5,4,4,4)+0.1);
 
 plot( out$s ~ out$t, col="blue", xlab="t(days)",ylab="percentage", ylim=c(0, 1.1), pch=16 );
 title ("rad52DD 050706 #1");                                        ############change here
 arrows( out$t, (out$s - out$e.s), out$t, (out$s + out$e.s), length=0.1, angle=90,code=3, lty=1, col="blue" );

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 lines( fit.s ~ t, col="blue",lty=2);
 
 points( out$Pb ~ out$t, pch=15 );
 arrows( out$t, (out$Pb - out$e.b), out$t, (out$Pb + out$e.b), length=0.1, angle=90,code=3, lty=1 );
 lines( fit.b ~ t, lty=2);

 ##calculate Rb = d(Pb)/dt
 # derivative.black <- function(b.max, b.min, T, w, t) {  (b.max - b.min) * w * t ^(w -1) / (1 + (t / T ))^2; }
 out$dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );  #first point is Inf
 out$g.e = genome.integrity(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );

 out$Rb     = out$dPb / out$g.e;                        ### rate of becoming blacks
 out$R0.5   = out$R0.5.raw / out$g.e;                   #### 090706 change
 ## lines ( out$Rb ~ out$t, col="brown", lty=2) 
 #lines ( out$R0.5 ~ out$t, col="green", lty=2); 
 #points ( out$R0.5 ~ out$t, col="green");

 # out$e.b05  =  ifelse( out$e.b05==0, out$R0.5 * (out$e.b/out$Pb), out$e.b05    );
 # out$e.dPb =  out$dPb * (out$e.b / out$Pb)  ## I use out$e.b as proxy  
 
 out$L   = out$R0.5 / out$Rb   #This is equivalent to out$R0.5.raw / out$dPb
 out$L[3] = 0;
 out$e.L = out$L * sqrt( (out$e.b05 / out$R0.5)^2 + (out$e.b / out$Pb)^2 )
 ylim = c(0, max(out$L[2:length(out$L)])*1.7 )

 par(new=T);
 plot( out$L ~ out$t, col="red", axes=F, xlab="", ylab="", ylim=ylim, pch=16 ); 
 arrows( out$t, (out$L - out$e.L), out$t, (out$L + out$e.L), length=0.1, angle=90,code=3, lty=1, col="red" );
 axis(4, at=pretty(ylim) )
 mtext( "ratio", 4, 2);

 L = out$L[2:length(out$L)]; t=out$t[2:length(out$L)];
 L = c( 0, L);  t = c(0, t);
 fm.L = gnls( L ~ A * exp(B*t) + C, start=list( A=0.005, B=0.5, C=0.005) );
 t = seq(0, max(out$t),by=0.1);
 fit.L = fm.L$coefficients[1] * exp( fm.L$coefficients[2] * t) + fm.L$coefficients[3]
 lines( fit.L ~ t, col="red", lty=2);

 #add labels here
 cols  =c("blue","black","red");
 labels=c("viability","Pr(b)","half/full")
 ltypes=c(2,2,2) 
 pch   =c( 16, 15,16)
 legend( 5 , 10 ,labels, col=cols, lty=ltypes, pch=pch);

# dev.off();

# quit("yes");


