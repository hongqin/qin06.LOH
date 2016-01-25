# only calculate T.c, T.g
rm(list=ls());

library(nlme)
file = "050306.M34.1.tab";

logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) {
   (b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w ); 
}

genome.integrity <- function(b.max, b.min, T, w, t) {
  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max);
}


tb= read.table( file, header=T, sep="\t", fill=T);
tb$t = tb$t / 24;
tb2 = tb;
labels = names( tb );

#normalize tb2
 tmp.indice = seq( 1 : length(labels) );
 names( tmp.indice ) = labels;
 for( j in tmp.indice["white"] : length(labels) ) {
 #for( j in 5:13 ) {
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
 col.labels = c( labels[ c(1,tmp.indice["white"] : length(labels)) ], "total" ); ######## bug here 091806
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
 par(new=T)
 plot( fit.m ~ t, col="brown", axe=F);
 
 #estimate errors for s
 error.s = out$s -  logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], out$t );
 error.s = abs(error.s)
 out$e.s = ifelse( out$e.s==0, error.s, out$e.s );
 out$e.s= ifelse( out$e.s<0.01, 0.01, out$e.s );
 
######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)
 b = out$Pb
 t = out$t;
 b.min = min( b ); b.min  # 0.008868349
 b.max = max( b ); b.max  # 0.055
 fb = function( t, T.g, w ) { b.max - (b.max - b.min) / (1 + (t/T.g)^w) }

# ws = 1 / out$e.b ^ 2
 #####   Pr[B>=b] = b.max - (b.max - b.min) /( 1 + ( t / T )^ w ) ########formula for cumulative Pr(b)
 #fm.b = gnls( b ~ 0.0633694 - (0.0633694 - 0.008868349) / (1 + (t/T.g)^w ), start=list( T.g=5, w=3), weights = ws);
 fm.b = gnls( b ~ fb(t, T.g, w), start=list( T.g=5, w=5));
 fm.b

 #estimate error for black
 error.b = out$Pb - logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], out$t);
 error.b = abs( error.b );
 out$e.b = ifelse( out$e.b==0, error.b, out$e.b );
 #out$e.b = ifelse( out$e.b<0.01, 0.01, out$e.b);

 plot( fm.b );

 t = seq(0, max(out$t),by=0.1);
 fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 
 plot( out$Pb ~ out$t );
 lines( fit.b ~ t, col="blue");

###################overlay s, Pb, 
 postscript("091706.100505M34.s.Pb.ps", width=8, height=8)
 par(mar=c(5,4,4,4)+0.1);
 
 plot( out$s ~ out$t, col="blue", xlab="t(days)",ylab="Viability", ylim=c(0, 1.1), pch=16 );
 title (file);                                        ############change here
 arrows( out$t, (out$s - out$e.s), out$t, (out$s + out$e.s), length=0.1, angle=90,code=3, lty=1, col="blue" );

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 lines( fit.s ~ t, col="blue",lty=2);
 
 par(new=T)
 fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 ylim = c( 0, 0.1) 
 plot( out$Pb ~ out$t, pch=15, axes=F,xlab="",ylab="", ylim=ylim );
 arrows( out$t, (out$Pb - out$e.b), out$t, (out$Pb + out$e.b), length=0.1, angle=90,code=3, lty=1 );
 lines( fit.b ~ t, lty=2);
 axis(4, at=pretty( ylim) )
 mtext( "P(b)", 4, 2);

 cols  =c("blue","black" );
 labels=c("viability","P(b)")
 ltypes=c(2,2,) 
 pch   =c(16,15)
 legend( 9.5 , 1.1 ,labels, col=cols, lty=ltypes, pch=pch);

 dev.off();

###################overlay s, m, Pb, Rb
 postscript("091706.100505M34.m.Rb.ps", width=8, height=8)
 par(mar=c(5,4,4,4)+0.1);
 
 #plot( out$s ~ out$t, col="blue", xlab="t(days)",ylab="percentage", ylim=c(-0.05, 1.2), pch=16 );
 #arrows( out$t, (out$s - out$e.s), out$t, (out$s + out$e.s), length=0.1, angle=90,code=3, lty=1, col="blue" );

 t = seq(0, max(out$t),by=0.1);

 fit.dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], t );  #first point is Inf
 fit.g   = genome.integrity(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], t );
 fit.Rb  = fit.dPb / fit.g;                        ### rate of becoming blacks
 plot( fit.Rb ~ t, col="black", xlab = "days", ylab="Black rate", type='l');
 title (file);

 par(new=T)
 fit.m = logistical.mortality ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.m ~ t, col="red", axe=F, xlab="",ylab="", type='l');
 axis(4, at=pretty(fit.m) )
 mtext( "Mortality rate", 4, 2);

 #cols  =c( "blue","black", "purple", "red" );
 #labels=c("viability","P(b)", "mortality rate", "black rate")
 cols   =c( "black", "red", "blue")
 labels =c("black rate", "mortality rate", "viability")
 ltypes =c(1,1,2,1) 
 pch    =c( NA, NA,NA,16,15,)
 legend( 5 , max(fit.m)/3, labels, col=cols, lty=ltypes, pch=pch);

 par(new=T);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue",lty=2, axes=F, xlab="",ylab="", type='l');

 #par(new=T)
 #fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 #points( out$Pb ~ out$t, pch=15 );
 #arrows( out$t, (out$Pb - out$e.b), out$t, (out$Pb + out$e.b), length=0.1, angle=90,code=3, lty=1 );
 #lines( fit.b ~ t, lty=2);
 #lines( fit.g ~ t, col="green")

 dev.off();


##################################
 postscript("091706.100505M34.Rb.R05.L.ps", width=8, height=8)
 par(mar=c(5,4,4,4)+0.1);

 ##calculate Rb = d(Pb)/dt
 # derivative.black <- function(b.max, b.min, T, w, t) {  (b.max - b.min) * w * t ^(w -1) / (1 + (t / T ))^2; }
 out$dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );  #first point is Inf
 out$g.e = genome.integrity(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );

 out$Rb     = out$dPb / out$g.e;                        ### rate of becoming blacks
 out$R0.5   = out$R0.5.raw / out$g.e;                   #### 090706 change
 ylim = c(0, 0.05)
 plot( out$R0.5 ~ out$t, col="green", xlab="days",ylab="half-black rate", type='l',lty=2, ylim=ylim  );
 title (file);                                        ############change here
 points( out$R0.5 ~ out$t, col="green", pch=15);
 arrows( out$t, (out$R0.5 - out$e.b05), out$t, (out$R0.5 + out$e.b05), length=0.1, angle=90,code=3, lty=1, col="green" );

 par(new=T);
 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue", type='l',lty=2, axes=F, xlab="",ylab=""  );

 #add labels here
 cols  =c("blue","black","red", "green");
 labels=c("viability","Rb","half/full", "R0.5")
 ltypes=c(2,2,1,1) 
 pch   =c( NA,NA,NA, 16)
 legend( 0 , 0.8 ,labels, col=cols, lty=ltypes, pch=pch);
 

 par(new=T);
 fit.dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], t );  #first point is Inf
 fit.g   = genome.integrity(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], t );
 fit.Rb  = fit.dPb / fit.g;                        ### rate of becoming blacks
 plot( fit.Rb ~ t, col="black", xlab ="", ylab="", type='l', axes=F, lty=2);
 axis(4, at=pretty(fit.Rb) )
 mtext( "black rate", 4, 2);

 par(new=T);
 out$L   = out$R0.5 / out$Rb   #This is equivalent to out$R0.5.raw / out$dPb
 out$L[1:3] = NA;
 #out$e.L = out$L * sqrt( (out$e.b05 / out$R0.5)^2 + (out$e.b / out$Pb)^2 )
 ylim = c(0, max(out$L[2:length(out$L)])*1.7 )

 plot( out$L ~ out$t, col="red", axes=F, xlab="", ylab="",  pch=16, type='l' ); 
 #arrows( out$t, (out$L - out$e.L), out$t, (out$L + out$e.L), length=0.1, angle=90,code=3, lty=1, col="red" );

 dev.off();

quit("yes");




 L = out$L[2:length(out$L)]; t=out$t[2:length(out$L)];
 L = c( 0, L);  t = c(0, t);
 fm.L = gnls( L ~ A * exp(B*t) + C, start=list( A=0.005, B=0.5, C=0.005) );
 t = seq(0, max(out$t),by=0.1);
 fit.L = fm.L$coefficients[1] * exp( fm.L$coefficients[2] * t) + fm.L$coefficients[3]
 lines( fit.L ~ t, col="red", lty=2);

 # out$e.b05  =  ifelse( out$e.b05==0, out$R0.5 * (out$e.b/out$Pb), out$e.b05    );
 # out$e.dPb =  out$dPb * (out$e.b / out$Pb)  ## I use out$e.b as proxy  
 


