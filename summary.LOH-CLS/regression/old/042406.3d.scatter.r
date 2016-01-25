 rm( list = ls() );
 tb = read.table( "data.tab", header=T);
 tb$L0 = tb$R0.5.0 / tb$Pblack.0;
 tb$Lmax = tb$R0.5.max / tb$Pb.max;

 attach(tb);
 library(scatterplot3d);

 postscript("L0.ARLS.CLS.3d.042406.ps", width=8,height=8); 
 #jpeg( "L0.ARLS.CLS.3d.042406.jpg");
 s3d <- scatterplot3d( ARLS, CLS, L0, type='h',color="red", pch=16,
  main="The correlations among ALRS, CLS, and L0", 
  );
 my.lm <- lm( tb$L0 ~ tb$ARLS + tb$CLS);
 s3d$plane3d( my.lm, lty.box="solid");
 text(s3d$xyz.convert( ARLS, CLS, L0), labels=tb$strain, pos=1)
 dev.off();

 postscript( "basal.L0.R0.5.Pblack.3d.ps", width=8,height=8 );
 s3d2 <- scatterplot3d( R0.5.0, Pblack.0, L0, type="h",color="blue",pch=16)
 m3 = lm( tb$L0 ~ tb$R0.5.0 + tb$Pblack.0 );
 s3d2$plane3d( m3 );
 text(s3d2$xyz.convert (R0.5.0, Pblack.0, L0), labels=tb$strain, pos=1);
 dev.off();
