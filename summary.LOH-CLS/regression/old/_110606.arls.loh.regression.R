 rm( list = ls() );
 tb = read.table("_110606.summary.loh.arls.csv", header=T, sep="\t");
 tb.old = tb;
 tb = tb[c(1:11), c(1:5,8, 10, 12, 14, 16, 18, 20, 22)]
 #tb[2:length(tb[1,]), ] = as.numeric( tb[2:length(tb[1,]), ] );

 summary( lm( ARLS ~ Tg, data=tb) ); #p = 0.36
 summary( lm( ARLS ~ Tc, data=tb) ); # p 0.48
 summary( lm( ARLS ~ Tmmax, data=tb) ); #p 0.59
 summary( lm( ARLS ~ Tbmax, data=tb) ); #p 0.36
 summary( lm( ARLS ~ Td, data=tb) ); #p 0.93

### begin ####ARLS ~ b.max 
 summary( lm( ARLS ~ b.max, data=tb) ); #p 0.06455
 
 long.arls  = tb$b.max[ tb$ARLS > 33  ]
 short.arls = tb$b.max[ tb$ARLS < 33  ]
 t.test( long.arls, short.arls)   # 0.02863
 
 tb.tmp = tb[, c("ARLS", "b.max")] 
 tb.tmp$by.arls = ifelse( tb$ARLS  >33, 1, 0)
 tb.tmp$by.b    = ifelse( tb$b.max > 0.16, 1, 0)

 # wrong:: fisher.test( tb.tmp[,c("by.arls","by.b")] )
 tmp = table( tb.tmp[,c("by.arls","by.b")] )
 fisher.test( tmp ); #p 0.06

 postscript("110606.arls.bmax.ps", width=8, height=8, horizontal=F)
 m1 = lm( b.max ~ ARLS, data=tb);
 plot( b.max ~ ARLS, data=tb, pch=16, col="black", main="b.max ~ ARLS", xlim=c(22,40), ylim=c(0.02, 0.3) );
 abline( m1 , col= "red");
 y = tb$b.max +0.01 ;
 names(y) = as.character( tb$strain );
 y[c("M32" )] = y[c("M32" )] - 0.02 
 x = tb$ARLS +0.5 ;
 x[c("M1-2")] = x[c("M1-2")] - 2
 text( x, y, tb$strain);
 dev.off();

###end #######ARLS ~ b.max 

### begin ALRS ~ (Tg-Tc)/Tc
 tb$frac.gc = (tb$Tg - tb$Tc) / tb$Tc;
 tb$frac.bm = (tb$Tbmax - tb$Tmmax) / tb$Tmmax; 
 
 summary( lm( ARLS ~ Tg:Tc, data=tb) );
 summary( lm( ARLS ~ frac.gc, data=tb)) # p 0.007
 summary( lm( ARLS ~ frac.bm, data=tb)) # p 0.01

 postscript("110606.Tg.Tc.arls.ps",width=8,height=8,horizontal=F);
 m2 = lm( frac.gc ~ ARLS, data=tb);
 plot( frac.gc ~ ARLS, data=tb, pch=16, col="black", main="(Tg-Tc/Tc) ~ ARLS", xlim=c(22,40) );
 abline( m2 , col= "red");
 y = tb$frac.gc + 0.01
 names(y) = as.character( tb$strain );
 #y[c("M32" )] = y[c("M32" )] - 0.02 
 x = tb$ARLS ;
 #x[c("M1-2")] = x[c("M1-2")] - 2
 text( x, y, tb$strain, pos=4);
 dev.off();


 postscript("110606.Tmax.m.g.arls.ps",width=8,height=8,horizontal=F);
 m3 = lm( frac.bm ~ ARLS, data=tb);
 plot( frac.bm ~ ARLS, data=tb, pch=16, col="black", main="(Tmax.g-Tmax.m/Tmax.m) ~ ARLS", xlim=c(22,40) );
 abline( m3 , col= "red");
 y = tb$frac.bm + 0.01
 names(y) = as.character( tb$strain );
 x = tb$ARLS ;
 text( x, y, tb$strain, pos=4);
 dev.off();
### end 






#quit("yes");
