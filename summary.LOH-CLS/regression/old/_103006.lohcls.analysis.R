 rm( list = ls() );
 tb = read.table( "loh.csv", header=T, sep="\t");
 tb = tb[1:37,]

 tb$L0 = tb$Half.0.raw / tb$B0;

 m1 = lm( Tc ~ Tg + expt + strain + Bmin + Bmax + B0 + Half.0.raw + L0, data=tb)
 summary( m1 ); #nonsense

 m2 = lm( Tc ~ Tg + strain + Bmin + Bmax + B0 + Half.0.raw + L0, data=tb)
 summary( m2 );

 m3 = lm( Tc ~ Tg + Bmin + Bmax + B0 + Half.0.raw + L0, data=tb)
 summary( m3 );

 m4 = lm( Tc ~ Tg , data=tb)
 summary( m4 );

quit("yes");
