rm(list=ls());

tb = read.table( "103006.lohcls.summary.csv", header=T, sep="\t");

#tb2 = t(tb);

strains = sort( unique( tb$strain) )

out.Tc = data.frame( matrix(  nrow=max(table(tb$strain)), ncol=length(strains) ) );
names(out.Tc) = strains;

out.Tg = out.Tc;

for( ss in 1:length(strains) ) {
  my.Tc = tb$Tc[tb$strain == strains[ss]];
  my.Tg = tb$Tg[tb$strain == strains[ss]];
  out.Tc[ 1:length(my.Tc),ss] = my.Tc;
  out.Tg[ 1:length(my.Tg),ss] = my.Tg;
}

x = as.vector( table(tb$strain) )

out2 = t ( rbind( as.character(strains), mean(out.Tc, na.rm=T), sd(out.Tc, na.rm=T), mean(out.Tg, na.rm=T), sd(out.Tg, na.rm=T) , x ) )


out2 = data.frame( out2 );

names(out2) = c("strain", "Tc", "Tc.sd", "Tg", "Tg.sd", "Nloh");

write.table( out2, "103106.Tc.Tg.summary.csv",  quote=F, sep="\t", row.names=F);

write.table( out.Tc, "_103106.Tc.values.csv", quote=F, sep="\t", row.names=F);

