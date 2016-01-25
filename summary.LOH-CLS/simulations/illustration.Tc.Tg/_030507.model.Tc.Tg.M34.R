### 030507 model Tc Tg for publication
### I want to explain Tc, Tg, and (Tg-Tc)/Tc

rm(list=ls());
library(nlme)

Tc = 6;
wc = 10;
wg = 10;
Tg = 8
b.min = 0.02
b.max = 0.25

######## functions
logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) {(b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w );}

genome.integrity <- function(b.max, b.min, T, w, t) {  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max); }
######### end of functions


 t = seq( 0, 15, by=0.05);
 header = c("t","Pb","s", "g","m","Rb","R0.5", "L"  );   ### L = rate(1/2) / rate(black) 

 out = data.frame( matrix( nrow= length(t) , ncol= length(header) ) );
 names( out ) = header;

 out$t = t;     # "t"

 #######calculate s , g
 out$s      = logistical.viability( Tc, wc, t);
 out$Pb     = logistical.black( b.max, b.min, Tg, wg, t);  
 out$g      = 1 - 2 * out$Pb;

 pdf( paste("030507", "modelTgTc", "pdf", sep="."), width=8, height=8 ); 

 plot( out$s ~ out$t , type='l', main= "", col="blue", xlab="Time (t)", ylab="Percentage",
      xlim=c(0,18)	)
 box();
 #axis(1, at = pretty(c(0,18)) );
 #axis(2, at= c(seq(0,1,by=0.2),0.5) );

 lines( out$g ~ out$t, col="red");

 labels = c("Viability (s)", "Genome integrity (g)");
 ltypes = c(1,1, 1);
 legend( (max(out$t)*0.7 ), 0.95, labels, col=c("blue", "red" ), lty = ltypes);

 points( Tc, 0.5, pch=19, col="blue", cex=1.2 )
 points ( Tg,  (1-b.max - b.min), pch=19, col="red", cex=1.2);

 arrows( Tc, 0.5, Tc, -1, lty=2, col="black" );
 arrows( Tg, (1-b.max - b.min), Tg, -1, lty=2, col="black" );
 mtext( "Tc",side=1,at=c(Tc) );
 mtext( "Tg",side=1,at=c(Tg) );

 arrows( Tc, 0.5, -5, 0.5, lty=2 );
 # text( -0.2, 0.47, "0.5" );
 arrows( Tg, (1-b.max - b.min), -5, (1-b.max - b.min), lty=2 );
 # mtext( (1-b.max - b.min), side=2, at=c((1-b.max - b.min)) );
 # text( 1, (1-b.max - b.min)*0.95, "1 - b.max - b.min" );


 dev.off();
