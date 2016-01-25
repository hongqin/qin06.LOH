 rm( list = ls() );
 tb = read.table( "loh.csv", header=T, sep="\t");
 tb$L0 = tb$half.0.raw / tb$b0;

 tb.old = tb;

 tb[as.character(tb$expt) == "M34.100506.37C.tab",    c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128.121205.37C.tab", c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128,053006",         c(1,3:length(tb[1,])) ] = NA;

 tb = tb[! is.na(tb$Strain), ]
 strains = sort( unique( tb$Strain) )

##to do here

#quit("yes");
