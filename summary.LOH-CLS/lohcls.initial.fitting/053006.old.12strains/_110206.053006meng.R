# 110606 add Td for half black increasing rate

rm(list=ls());
library(nlme)

file = "053006.clsloh.old.12strains.tab";

tb0= read.table( file, header=T, sep="\t", fill=T);

strains = as.character( unique( tb0$Strain ) );

######## functions
logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) {(b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w );}

genome.integrity <- function(b.max, b.min, T, w, t) {  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max); }
############ end of functions

  # "101S""M1-2" M13 M14" M2-8-1" "M2-8-2" "M32"    "M34" "M5"  "M8" "YPS128" "YPS163"
ww = c(12,17,17,       5,   10,    10,       6,      10,   20,  10,     5,   5)
TT = c(4,  7,   10,   7,    6,        4,       6,     5,   5,   10,     4,   4);

parameters = c("strain", "Tc","wc", "Tmmax","Tg", "wg", "Tbmax","Td", "wd", "Tdmax","b.min", "b.max", "b0", "half.0.raw");
out2 = data.frame( matrix( nrow= length(strains), ncol=length(parameters) ) );
names(out2) = parameters;

pdf( paste("110206", file, "pdf", sep="."), width=8, height=8 ); 

###################################
for( ss in 1:length(strains) ) {
# ss =1; 
 out2$strain[ss] =  strains[ss]; 
 current_strain =   strains[ss];

 tb = tb0[as.character(tb0$Strain)==current_strain, 2:length(tb0[1,])];
 #if(  current_strain == "M5" ) { tb[c(4), ] = NA;   }
 tb[ is.na(tb$white), ] = NA;
 tb = tb[ ! is.na(tb$t),  ]

 tb$t = tb$t / 24;

 labels = names(tb);
 tmp.indice = seq( 1 : length(labels) );
 names( tmp.indice ) = labels;

 for( i in 1:length(tb[,1])) { #row
   for( j in tmp.indice["white"] : length(labels) ) {
      if (is.na( tb[i,j])) { tb[i,j] = 0; }
   }
 }

 #normalization
 tb2 = tb;  
 for( j in tmp.indice["white"] : length(labels) ) {
   tb2[,j] = tb2[,j] * tb2$D1 * tb2$D2;
 }

 tb.sub = tb2 #for minimal coding 

 #generate row indice for averaging
 row.steps = as.vector( table( tb.sub$t ) ); # a new trick, ha

 #generate tb.m
 n.row = length( row.steps );
 col.labels = c("t", "white", "black", "B0.5","B0.25", "B0.125","B0.0625","B2", "Btiny","total" );

 tb.m = data.frame( matrix( nrow=n.row, ncol= length(col.labels) ) )  # mean values
 names( tb.m ) = col.labels; 
 
 upper.row = 0; #up pointer
 lower.row = 0; #low pointer 
 for( i in 1:n.row ) {
   upper.row = lower.row + 1;
   lower.row = upper.row + row.steps[ i ] - 1 ;

   tb.m[i,1]  = tb.sub$t[ upper.row ]
   for( j in 2: ( length(col.labels) - 1 ) ) {
     tb.m[i, j] = mean( tb.sub[ upper.row : lower.row, j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:( length(col.labels) - 1 )], na.rm=T );
 }

 header = c("t","Pb","s", "g","m","Rb","R0.5", "L"  );   ### L = rate(1/2) / rate(black) 
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

 #######calculate s , g
 out$s      = tb.m$total / tb.m$total[1]
 out$Pb     = tb.m$black / tb.m$total;
 out$g      = 1 - 2 * out$Pb;
  # out$e.s    =  out$s * (tb.m$sd.w / tb.m$white);
  # out$e.b    =  out$Pb * (tb.m$sd.b / tb.m$black);

 out$R0.5.raw   = tb.m$B0.5 /tb.m$total; #needs to be adjusted by g.e
  # out$e.b05  =  out$R0.5.raw * ( tb.m$sd.b05 / tb.m$B0.5)

 ### plot of the raw data 
 plot( out$s ~ out$t , type='l', main= paste(current_strain), col="blue", log='y');
 lines( out$Pb ~ out$t, col="black");
 labels = c("viability","black");
 ltypes = c(1,1);
 legend( (max(out$t)*0.7 ), 0.8, labels, col=c("blue", "black"), lty = ltypes);

 ### estimate half life
 t1 = max( out$t[out$s>0.5] )
 t2 = min( out$t[out$s<0.5] )
 t.start = (t1+t2)/2 

 ######## calclulate T0.5 w_s
 t= out$t
 s= out$s
 fs2 = function(t,v) {  1/( 1 + ( t / v )^ ww[ss] ) }

 #if ( strains[ss] == "M13") { s[c(2,3,4)] = 1 }

 if (( strains[ss] == "M14") | (strains[ss]=="M13" ) | (strains[ss]=="M2-8-1" )) {
  fm.s <- gnls( s ~ fs2(t,v) , start = list( v = t.start)  );
  out2[ss, c("Tc", "wc")] = c( fm.s$coefficients[c("v")], ww[ss] );
 } else {
  fm.s <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = t.start, w = ww[ss] )  );
  fm.s # this the half life T1/2 w_s
  out2[ss, c("Tc", "wc")] = fm.s$coefficients[c("v","w")]
 }

 tt = seq(0, 2*max(out$t),by=0.05);
 fit.m = logistical.mortality ( out2$Tc[ss], out2$wc[ss], tt );
 out2$Tmmax[ss] = tt[fit.m == max(fit.m)]


######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)
 out$Pb = ifelse( out$Pb==0 & out$t > t2, NA, out$Pb);   #remove outliers of Pb in late stage

 if ( strains[ss] == "101S")  {  out$Pb[c(4)] = NA;   }
 if (strains[ss] =="M2-8-2") { out$Pb[c(5)] = 0.2; }

 #estimate starting values
 b.tmp = sort( out$Pb[ out$Pb>0 ] );
 b.min = b.tmp[1]/2 + b.tmp[2]/2;
 b.tmp = rev(b.tmp); 
 b.max = b.tmp[ 2 ] /2 + b.tmp[ 2] /2
 # b.min = b.tmp[1]/2 + b.tmp[2]/2;
 #  b.max = b.tmp[ length(b.tmp) ] /2 + b.tmp[ length(b.tmp) -1 ] /2

 b = out$Pb[ ! is.na(out$Pb) ]
 t = out$t[  ! is.na(out$Pb) ]  ;

 fb = function( t, T.g, w ) { b.max - (b.max - b.min) / (1 + (t/T.g)^w) }
 fb2 = function( t, T.g)  { b.max - (b.max - b.min) / (1 + (t/T.g) ^ww[ss] ); }

 #####   Pr[B>=b] = b.max - (b.max - b.min) /( 1 + ( t / T )^ w ) ########formula for cumulative Pr(b)
 #fm.b = gnls( b ~ b.max2 - (b.max2 - b.min2) / (1 + (t/T.g)^w ), start=list( T.g= t2, w=10) );
 if( ( strains[ss]=="M13")|(strains[ss]=="M2-8-2") |(strains[ss]=="M8") ) {
  fm.b = gnls( b ~ fb2( t, T.g), start=list( T.g= t2) );
  fm.b
  out2[ss, c("Tg", "wg")] = c( fm.b$coefficients[c("T.g")], ww[ss]);
 } else {
  fm.b = gnls( b ~ fb( t, T.g, w), start=list( T.g= t2, w= ww[ss]) );
  fm.b
  out2[ss, c("Tg", "wg")] = fm.b$coefficients[c("T.g","w")]
 }
 out2[ss, c("b.min", "b.max", "b0", "half.0.raw")] = c(b.tmp[1], b.tmp[length(b.tmp)], out$Pb[1], out$R0.5.raw[1]  );

 fit.dPb = derivative.black(b.max, b.min, out2$Tg[ss], out2$wg[ss], tt );  #first point is Inf
 fit.g   = genome.integrity(b.max, b.min, out2$Tg[ss], out2$wg[ss], tt );
 fit.Rb  = fit.dPb / fit.g;     

 out2$Tbmax[ss] = tt[fit.Rb == max(fit.Rb)]


 ######## calculate Td, for half blacks
 out$g.e = genome.integrity(b.max, b.min, out2$Tg[ss], out2$wg[ss], out$t );
 out$R0.5   = out$R0.5.raw / out$g.e;                   ##

 out$dPb = derivative.black(b.max, b.min, out2$Tg[ss], out2$wg[ss], out$t );  #first point is Inf
 out$Rb     = out$dPb / out$g.e;   

 out$L     = out$R0.5 / out$Rb;   

 plot( out$R0.5 ~ out$t, main= paste(current_strain, "050806"), type='b' );
 plot( out$Rb   ~ out$t, main= paste(current_strain, "050806"), type='b' );
 plot( out$L    ~ out$t, main= paste(current_strain, "050806"), type='b', log='y' );

 #if( e[ss]== "M5.032906.tab" ) { out$R0.5[7]=NA;  }
 #if( expts[ss]== "YPS128.121205.37C.tab" ) { out$R0.5[6]=NA;  }

 out$R0.5[ out$R0.5==0 & out$t > t1 ] = NA; 
 if (strains[ss] =="M2-8-1" ) { out$R0.5[c(4,6)] =c(0.01, out$R0.5[5] )  }

 half.tmp = sort( out$R0.5[ out$R0.5>0 ] );
 half.min = half.tmp[1]/2 + half.tmp[2]/2;  #the first may be an outlier
 half.tmp = rev(half.tmp); 
 half.max = half.tmp[ 1 ] /2 + half.tmp[ 2 ] /2

 half = out$R0.5[ ! is.na(out$R0.5) ]
 t    = out$t[  ! is.na(out$R0.5) ]  ;
 fd = function( t, T.d, w ) { half.max - (half.max - half.min) / (1 + (t/T.d)^w) }
 fd2 = function( t, T.d)  { half.max - (half.max - half.min) / (1 + (t/T.d) ^ww[ss] ); }

 if (strains[ss] == "M2-8-2") { 
   fm.h = NA;
   out2[ss, c("Td", "wd")] = c( NA, NA );
 } else if (strains[ss] == "101S") { 
   fm.h = gnls( half ~ fd2( t, T.d), start=list( T.d= TT[ss] ) );
   fm.h
   out2[ss, c("Td", "wd")] = c( fm.h$coefficients[c("T.d")], ww[ss] );
 }else {
   fm.h = gnls( half ~ fd( t, T.d, w), start=list( T.d= t2, w= ww[ss]) );
   fm.h
   out2[ss, c("Td", "wd")] = fm.h$coefficients[c("T.d","w")]
 }

 #use weighted average to estimate the peak of R0.5
 tmp = max(out$R0.5, na.rm=T); 
 td.max = out$t[ out$R0.5 == tmp ]
 td.max = td.max[! is.na(td.max)]
 out$tmp = seq(1, length(out$t));
 td.max.pos = out$tmp[ out$t == td.max ]
 max.range = c( td.max.pos -1, td.max.pos );
 if( td.max.pos < length(out$t) ) {
    max.range = c( max.range, td.max.pos +1);
 } else {
    max.range = c( td.max.pos -2, max.range);
 }
 y = out$R0.5[max.range];
 y[is.na(y)] = 0;
 out2$Tdmax[ss] = sum( out$t[max.range] * y ) / sum(y)  #??????????

 write.table(out2, "estimations.meng053006.110206.csv", row.names=F, quote=F, sep="\t");

}#end of ss loop

 write.table(out2, "estimations.meng053006.110206.csv", row.names=F, quote=F, sep="\t")

 dev.off();
