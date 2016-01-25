
rm(list=ls());

tb = read.table("time.zero.all.csv", header=T,fill=T,sep="\t")
tb.old = tb;

tb = tb.old[ ! is.na(tb$t), ]

strains = as.character( levels(tb$Strain) );
strains = strains[ 2: length(strains) ]; #the first is null

out = data.frame( strains);

for( ss in 1:length(strains) ) {
 sub = tb[ tb$Strain == strains[ss],  ]
 sub$black = ifelse( is.na(sub$black), 0, sub$black);  
 sub$B0.5  = ifelse( is.na(sub$B0.5),  0, sub$B0.5);  
 out$L0.all[ss]   = sum( sub$B0.5 ) / sum( sub$black );
 out$L0.small[ss] = sum( sub$B0.5[ !is.na(sub$white) ] ) / sum( sub$black[!is.na(sub$white)] );

 out$subtot.0.5[ss]   = sum( sub$B0.5 );  #013007 change
 out$subtot.black[ss] = sum( sub$black );  #013007 change

 my.tot = sum( sub[!is.na(sub$white), c("white", "black", "B0.5")] )
 out$Pb[ss] =   sum( sub$black[!is.na(sub$white)] ) / my.tot;
 out$Pb0.5[ss] = sum( sub$B0.5[!is.na(sub$white)] ) / my.tot;

 out$Pb.2[ss]  = mean( sub$black / sub$white, na.rm=T )
 out$Pb.sd[ss] = sd( sub$black / sub$white, na.rm=T )

 out$Pb0.5.2[ss]  = mean( sub$B0.5 / sub$white, na.rm=T )

 #my.half  = sub$B0.5[1];
 #my.white = sub$white[1]; 
 #for( j in 2:length(sub[,1]) ) {
 #  if( sub$B0.5[j] == 0 ) {
 #      my.white[length(my.half)] = my.white[length(my.half)] + sub$white[j] 
 #  } else {
 #      my.half  = c( my.half,  out$B0.5[j] );
 #     my.white = c( my.white, out$white[j] );
 #  }
 #}
 #out$Pb0.5.sd[ss] = sd( my.half / my.white);

 out$Pb0.5.sd[ss] = sd( sub$B0.5   / sub$white, na.rm=T ) #This put too much weight on zero half count

 #this does not work
 #out$Pb.sd[ss] = sqrt(  out$Pb[ss] * (1-out$Pb[ss]) );  
 #out$Pb0.5.sd[ss] = sqrt(  out$Pb0.5[ss] * (1-out$Pb0.5[ss]) ); 
}

 out$L0.sd = sqrt( (out$Pb0.5.sd / out$Pb0.5)^2 + (out$Pb.sd / out$Pb)^2 ) * out$L0.small;   

 #out$strains = as.character( out$strains);
 write.table( out, "013007.L0.Pb0.half.by.strain.csv", row.name=F, quote=F, sep="\t");



