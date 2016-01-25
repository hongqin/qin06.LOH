 rm( list = ls() );
 tb = read.table("summary.new.by.strain.csv", header=T, sep="\t");
 tb.old = tb;
 labels = names( tb.old );
 #tb = tb.old[c(1:11), c(1:5,8, 10, 12, 14, 16, 18, 20, 22)]

 tb = tb.old[c(1:11), c("strain","ARLS","R0","G","CLS","Tc", "Tg","Tmmax","Tbmax", "Td", "Tdmax","TLmax","Lmax", 
 "b.max", "b.min", "strains", "L0.all", "L0.small" , "Pbt0","Pb0.5t0", "Pbt0.b") ];


 tb$CLS.vs.Tc = tb$CLS / tb$Tc

 tb$Tg.vs.Tc = tb$Tg / tb$Tc

 save.image("021207.mean.by.strain.RData");

######### L0 ~ CLS
summary( lm( tb$L0.all ~ tb$CLS) ) #p 0.011  R2=0.53 <--

summary( lm( tb$L0.all ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G ) ) # none
# only CLS shows connection

summary( lm( tb$L0.small ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G ) ) #none

summary( lm( tb$L0.all ~ tb$ARLS) ) #p 0.49

summary( lm( tb$L0.all ~ tb$ARLS + tb$CLS) ) #p 0.05, L0~CLS partial is p=0.02


summary( lm( tb$L0.all ~ tb$Tc) ) #p 0.33

summary( lm( tb$L0.all ~ tb$Tg) ) #p 55           



####### b.max ~ ARLS
summary( lm( tb$b.max ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all ) ) #p 0.50

summary( lm( tb$b.max ~ tb$ARLS ))
#Residual standard error: 0.06868 on 9 degrees of freedom
#Multiple R-Squared: 0.4122,     Adjusted R-squared: 0.3469
#F-statistic: 6.311 on 1 and 9 DF,  p-value: 0.0332


######## Lmax
summary( lm( tb$Lmax ~ tb$ARLS ))  #p 0.98


summary( lm( tb$TLmax ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all + tb$b.max + tb$b.min + tb$Tdmax ) ) #none

summary( lm( tb$Lmax ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all + tb$b.max + tb$b.min + tb$Tdmax ) )
# 0.05 ???  but one degree of freedom?? So, this is questionable.


##### Tg/Tc ~ ARLS
summary( lm( Tg.vs.Tc ~ tb$ARLS) );
# Residual standard error: 0.1861 on 9 degrees of freedom
# Multiple R-Squared: 0.5551,     Adjusted R-squared: 0.5057
# F-statistic: 11.23 on 1 and 9 DF,  p-value: 0.00851



summary( lm( tb$CLS.vs.Tc ~ tb$Tg.vs.Tc ))  #p 0.35

summary( lm( tb$Tg.vs.Tc ~ tb$Lmax + tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all + tb$b.max + tb$b.min + tb$Tdmax ) )
## 0.05 ???


q("no");



postscript("110706.L0-CLS.ps", width=8,height=8, horizontal=F)
m.L0 = lm( tb$L0.all ~ tb$CLS)
plot(tb$L0.all ~ tb$CLS, xlab="CLS",ylab="L0",main="L0 ~ CLS", pch=16, xlim=c(2,18), ylim=c(0.05,0.26) )
abline( m.L0, col="red");

x = numeric( length(tb$CLS) ); names(x) = as.character( tb$strains ); 
y =x;
y[c("M34", "M5",  "M32", "M8")] = c(0.002, -0.01, -0.01, -0.01, )
strains = as.character(tb$strain); names(strains) = tb$strain
#pos = x; pos=NA;
strains["M5"] = "     M5";
text( tb$CLS, tb$L0.all+y+0.004, strains);

dev.off();

#########110606Mon regression:
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
