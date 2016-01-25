 #050707 change pdf size to 6x6
 #030707 
 rm( list = ls() );
 
 tb = loh;
 tb[tb$strain=="RAD52DD", 2] = NA;
 tb[tb$strain=="BY4743", 2] = NA;
 tb[tb$strain=="SGU57", 2] = NA;
 tb$expt[as.character(tb$expt) == "M34.100506.37C.tab"] = NA;
 tb$expt[as.character(tb$expt) == "YPS128.121205.37C.tab"] = NA;
# tb$expt[as.character(tb$expt) == "YPS128,053006"] = NA;

 tb2 = tb;
 tb = tb2[ ( ! is.na(tb$expt) ),]
 

 m.cg = lm( tb$Tg ~ tb$Tc )
 tmp1 = summary(m.cg);
 summary(m.cg); 
 pdf( "030707.Tc.Tg.pdf", width=6, height=6);
 plot( Tg ~ Tc, data=tb, pch=16
	# main="Tc~Tg, p=E-5, R^2=0.44"
  );
 abline( m.cg, col="red");
 text(8,16, "R^2=0.44,p<0.001")
 dev.off();

 m.bm = lm( tb$Tbmax ~ tb$Tmmax )
 summary(m.bm); #p 0.0008 R^2 = 0.30
 pdf( "030707.Tbmax.Tmmax.pdf", width=6, height=6 );
 plot( Tbmax ~ Tmmax, data=tb, pch=16 
	# main="Tbmax~Tmmax, p=1E05, R^2=0.44"
 );
 abline( m.bm, col="red");
 text(10, 15, "R^2=0.44, p<0.001")
 dev.off();






#quit("yes");
