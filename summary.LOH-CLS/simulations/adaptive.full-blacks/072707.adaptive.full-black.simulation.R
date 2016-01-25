### 072707 simulation to address adaptive full-black argument ###

rm(list=ls());
library(nlme)

#####key parameters; may change during simulations
 Nb0 = 1E5;
 Nw0 = Nb0 * 100; 
 Tc = 6;
 wc = 10;
 wg = wc;
 b.min = 0.01;
 b.max = 0.25;

 step = 0.05; 
 ##key parameters that may change during simulations
 Tg = NA;

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

#######Start, Model A: no dying of initial full blacks
 t = seq( 0, 15, by= step);
 header = c("t", "s", "Fw", "b", "Fb", "g","m","Rb","R0.5", "L"  );   

 out = data.frame( matrix( nrow= length(t) , ncol= length(header) ) );
 names( out ) = header;

 out$t = t;     # "t"
 out$s      = logistical.viability( Tc, wc, t);
 out$Fw = Nw0 * out$s;

 plot( out$s ~ out$t, type='l', col="blue", log='y', main="no dying of full-black" );
 
 #par( new = T );
 out$Fb = Nb0;
 out$b = out$Fb / (out$Fb + out$Fw);
 lines( out$b ~ out$t );

 text( 4, 1E-3, "no dying of full-black")
 text( 4, 1E-3/1.5, "leads to b->1");

 #When no black is dying, b(t) eventually will reach to 100%. This is obvious, in retrospect. 
#######End, Model A;

#######Start, Model B: dying of full-blacks is slower than white cells
 t = seq( 0, 15, by= step);
 header = c("t", "s", "Fw", "b", "Fb", "g","m","Rb","R0.5", "L"  );   

 out = data.frame( matrix( nrow= length(t) , ncol= length(header) ) );
 names( out ) = header;

 out$t = t;
 out$s = logistical.viability( Tc, wc, t);
 out$Fw = Nw0 * out$s;

 pdf ("072707.slower.dying.of.blacks.pdf", width=8, height=8);
 plot( out$s ~ out$t, type='l', col="black", log='', main="slower dying of full-blacks", xlim=c(0,22),
  xlab= "t", ylab="percentage" 
 );
 
 #par( new = T );
 out$Fb = Nb0;

 changes = c( 1, 1.25, 1.5, 2, 3, 10);
 colors = rainbow( length(changes) );
 for( i in 1:length(changes) ) {
   Tg = changes[i] * Tc; 
   tmp.b = logistical.viability( Tg, wc, t);  # slower dying of blacks;
   out$b = tmp.b * Nb0 / ( Nb0 + out$Fw);
   # out$g      = 1 - 2 * out$Pb;
   lines( out$b ~ out$t, col=colors[i] );
 }

 tmp.text = paste( "Tg=",changes,"Tc");
 labels = c("s", tmp.text );
 ltypes = c(1,1);

 legend( (max(out$t)*0.9 ), 0.7, labels, col=c("black", colors ), lty =ltypes );
 dev.off();
#######End, Model B;

q("yes");

 

