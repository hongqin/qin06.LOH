tb = read.table( "data.tab", header=T);
 tb$L0 = tb$R0.5.0 / tb$Pblack.0;
 tb$Lmax = tb$R0.5.max / tb$Pb.max;
 attach(tb);

 postscript("Pbmax-ARLS.042106.ps", width=8,height=8); 
 plot( Pb.max ~ ARLS, ylab="Maximal P(full-black)", xlim=c(24,38),ylim=c(0,0.35),
    main="ARLS is correlated with the maximal LOH"
 );
 xx = (ARLS+ 0.8); yy= Pb.max +0.01; yy[9]=yy[9]-0.02;xx[9]=xx[9]-0.5; yy[6]=yy[6]*0.9;
 text( xx, yy, tb$strain);

 my.lm <- lm( Pb.max ~ ARLS);
 abline( my.lm, col="red");

 dev.off();

