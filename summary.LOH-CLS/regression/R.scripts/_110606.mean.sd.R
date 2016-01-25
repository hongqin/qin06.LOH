 rm( list = ls() );
 tb = read.table( "loh.csv", header=T, sep="\t");
 tb$L0 = tb$half.0.raw / tb$b0;

 tb.old = tb;

 tb[as.character(tb$expt) == "M34.100506.37C.tab",    c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128.121205.37C.tab", c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128,053006",         c(1,3:length(tb[1,])) ] = NA;


 tb = tb[! is.na(tb$Strain), ]
 strains = sort( unique( tb$Strain) )

parameters =c("strain", "Nloh", "Tc", "Tc.sd", "Tg", "Tg.sd", "Tmmax","Tmmax.sd", "Tbmax","Tbmax.sd", 
"Td", "Td.sd", "Tdmax", "Tdmax.sd", "b.max", "b.max.sd", "b.min", "b.min.sd" );
out= data.frame( matrix(  nrow=length(strains), ncol=length(parameters) ) );
names(out) = parameters;
row.names(out) = strains;

out$strain = strains;

for( ss in 1:length(strains) ) {
  my.Tc = tb$Tc[tb$Strain == strains[ss]];
  out$Tc[ss]     = mean( my.Tc );
  out$Tc.sd[ss]  =   sd( my.Tc );
  out$Nloh[ss] = length( my.Tc); 

  my.Tg = tb$Tg[tb$Strain == strains[ss]];
  out$Tg[ss]     = mean( my.Tg );
  out$Tg.sd[ss]  =   sd( my.Tg );

  my.Tmmax = tb$Tmmax[tb$Strain == strains[ss]];
  out$Tmmax[ss]    =   mean( my.Tmmax);
  out$Tmmax.sd[ss] =     sd( my.Tmmax);

  my.b.max = tb$b.max[tb$Strain == strains[ss]];
  out$b.max[ss]    =   mean( my.b.max);
  out$b.max.sd[ss] =     sd( my.b.max);

  my.b.min = tb$b.min[tb$Strain == strains[ss]];
  out$b.min[ss]    =   mean( my.b.min);
  out$b.min.sd[ss] =     sd( my.b.min);

  my.Tbmax = tb$Tbmax[tb$Strain == strains[ss]];
  out$Tbmax[ss]    =   mean( my.Tbmax);
  out$Tbmax.sd[ss] =     sd( my.Tbmax);

  if( ( strains[ss] != "Rad52DD" )&( strains[ss] != "SGU57") ) {
  my.Td = tb$Td[tb$Strain == strains[ss]];
  out$Td[ss]     = mean( my.Td );
  out$Td.sd[ss]  =   sd( my.Td );

  my.Tdmax = tb$Td[tb$Strain == strains[ss]];
  out$Tdmax[ss]    =   mean( my.Tbmax);
  out$Tdmax.sd[ss] =     sd( my.Tbmax);
  }


}


write.table( out, "110606.loh.mean.sd.csv",  quote=F, sep="\t", row.names=F);


#quit("yes");
