### study change of mitotic asymmetry ###
rm(list=ls());
library(nlme)

######## functions
logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) 
{ ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) 
{(b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w );}

genome.integrity <- function(b.max, b.min, T, w, t) 
{  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max); }
####### end of functions

#####key parameters; may change during simulations
 Nb0 = 1E5;
 Nw0 = Nb0 * 100; 
 Tc = 3;
 wc = 2;
 wg = wc;

 step = 0.05; 
 ##key parameters that may change during simulations


#######Start, Model A, bfull and b05 share same Tg, 
 kfull = 5; kd = 100; 
 b.max = 0.4;
 b.min = b.max / kfull;
 d.max = 0.1;
 d.min = d.max / kd;
 Tg = 5;

 t = seq( 0, 15, by= step);
 header = c("t", "s", "bfull", "b05",  "L"  );   

 out = data.frame( matrix( nrow= length(t) , ncol= length(header) ) );
 names( out ) = header;

 out$t = t;     # "t"
 out$s = logistical.viability( Tc, wc, out$t );

 plot( out$s ~ out$t, type='l', col="blue", main="model A", ylim=c(0,1) );
 
 out$bfull = logistical.black(b.max, b.min, Tg, wg, out$t ) 
 lines( out$bfull ~ out$t );

 out$b05   = logistical.black(d.max, d.min, Tg, 10*wg, out$t ) 
 lines( out$b05 ~ out$t, col="green");

 out$L = out$b05 / out$bfull;
 lines( out$L ~ out$t, col="red");


#######End, Model A;


q("yes");

 

