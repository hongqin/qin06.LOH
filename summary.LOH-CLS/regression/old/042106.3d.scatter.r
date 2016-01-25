tb = read.table( "data.tab", header=T);
 tb$L0 = tb$R0.5.0 / tb$Pblack.0;
 tb$Lmax = tb$R0.5.max / tb$Pb.max;
 attach(tb);

 library(scatterplot3d);

 postscript("3D-scatter.042106.ps", width=8,height=8); 
 s3d <- scatterplot3d( ARLS, CLS, L0, type='h',color="red", pch=16,
  main="The correlations among ALRS, CLS, and L0", 
  );
 my.lm <- lm( tb$L0 ~ tb$ARLS + tb$CLS);
 s3d$plane3d( my.lm, lty.box="solid");
 dev.off();

