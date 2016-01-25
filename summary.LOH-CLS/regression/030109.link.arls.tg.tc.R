

 pdf("_030109.Tc.Tg.arls.pdf",width=6,height=6);
 plot( 1/tb$Tg.vs.Tc ~ tb$ARLS, pch=16, col="black", xlim=c(23,40), ylim=c(0.5,1.2),xlab="ARLS",ylab="Tc/Tg" );
 m2 = lm( 1/tb$Tg.vs.Tc ~ tb$ARLS );
 abline( m2 , col= "black", lty=2);
 y = 1/tb$Tg.vs.Tc + 0.01
 names(y) = as.character( tb$strain );
 x = tb$ARLS ;
 text( x, y, tb$strain, pos=4);
 text( 35, 1.15, "R2=0.52, p=0.012")
 dev.off();


 library(scatterplot3d);

 pdf("_030109.3D.ARLS.Tc.Tg.pdf", width=6,height=6); 
 ARLS = tb$ARLS
 Tc = tb$Tc
 Tg = tb$Tg
 s3d <- scatterplot3d( ARLS, Tg, Tc, type='h',color="red", pch=16,
  main="The correlations among ALRS, Tc, and Tg", angle=40, 
  #xlim=c(0,15),ylim=c(0,15),
  #zlim=c(25,40)
  );
 my.lm <- lm( Tc ~ ARLS + Tg  );
 s3d$plane3d( my.lm, lty.box="solid", col="blue" );
 dev.off();

 pdf("_030109.ARLS.Tc.partial.pdf", width=6,height=6); 
 ARLS = tb$ARLS
 Tc = tb$Tc
 Tg = tb$Tg
 m = lm( Tc ~ Tg );
 Tc2 = m$residuals

 plot( Tc2 ~ ARLS, pch=16, col="black", xlim=c(23,40),xlab="ARLS",ylab="Tc|Tg" );
 m2 = lm( Tc2 ~ ARLS );
 abline( m2 , col= "black", lty=2);
 y = Tc2 + 0.01
 x = ARLS ;
 y[9] = y[9] + 0.1; 
 y[11] = y[11]-0.15; x[11] = x[11]-2; 
 names(y) = as.character( tb$strain );
 text( x, y, tb$strain, pos=4)
 text( 35, 1.15, "R2=0.35, p=0.054")
 dev.off();

#q("no");

