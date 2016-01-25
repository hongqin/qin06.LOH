 rm( list = ls() );
 tb = read.table( "loh.csv", header=T, sep="\t");
 tb$L0 = tb$half.0.raw / tb$b0;

 tb.old = tb;
 tb = tb.old [1:37, ]

 tb[as.character(tb$expt) == "M34.100506.37C.tab",    c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128.121205.37C.tab", c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128,053006",         c(1,3:length(tb[1,])) ] = NA;

 tb = tb[! is.na(tb$Strain), ]

 m.cg = lm( tb$Tg ~ tb$Tc )
 summary(m.cg); #p -.0015 R^2 = 0.27
 postscript( "110706.Tc.Tg.ps", width=8, height=8, horizontal=F);
 plot( Tg ~ Tc, data=tb, pch=16,main="Tc~Tg, p=0.0015, R^2=0.27");
 abline( m.cg, col="red");
 dev.off();

 m.bm = lm( tb$Tbmax ~ tb$Tmmax )
 summary(m.bm); #p 0.0008 R^2 = 0.30
 postscript( "110706.Tbmax.Tmmax.ps", width=8, height=8, horizontal=F);
 plot( Tbmax ~ Tmmax, data=tb, pch=16,main="Tbmax~Tmmax, p=0.0008, R^2=0.30");
 abline( m.bm, col="red");
 dev.off();






#quit("yes");
