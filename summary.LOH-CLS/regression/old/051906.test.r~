 rm( list=ls() );
 library(scatterplot3d);

 tb = read.table( "data.tab", header=T);
 tb$L0 = tb$R0.5.0 / tb$Pblack.0;
 tb$Lmax = tb$R0.5.max / tb$Pb.max;
 
 
 my.lm <- lm( tb$L0 ~ tb$ARLS + tb$CLS);
 
 attach(tb);

 s3d <- scatterplot3d( ARLS, CLS, L0, type='h',color="red", pch=16,
  main="The correlations among ALRS, CLS, and L0",
  );
 s3d$plane3d( my.lm, lty.box="solid");

 text(s3d$xyz.convert (ARLS, CLS, L0), labels=tb$strain, pos=1);


 tb2 = tb[,c(1,2,5,10,12)]
 m <- lm( tb2$L0.2 ~ tb2$ARLS + tb2$CLS); summary(m);

