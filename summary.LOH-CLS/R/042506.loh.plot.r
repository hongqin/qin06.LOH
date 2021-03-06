# output basal loh and max loh

rm(list=ls());
postscript("042506.loh.plot.12strains.ps",width=8, height=8);

filelist = "_normalized.data.file";
tb = read.table(filelist, sep="\t");
for( kk in 1:length(tb[,1]) ) {  ###the file loop ###
 # kk=1;
 file = paste( "normalized.data/", as.character(tb[kk,1]), sep="" );
 #file = "INPUTFILE";
 #file = "normalized.data/norm.M28.022006.txt";

tb.m = read.table( file, header=T, sep="\t");

# todo: remove NA and zeros from tb.m
# output to out3,
# columns in out3:
 header = c("t","half.over.black","Pb","R0.5","s" );
 out425 = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out425 ) = header;

 out425[,1] = tb.m[,1]; # "t"

 #calculate the ratios and s
  for( i in length(out425[,1]):1 ) { #row
    out425$s[i] = tb.m$total[i]    / tb.m$total[1]; #s
    out425[i,2] = tb.m$B0.5[i]     / tb.m$black[i];   
  }

 #calculate the risks (in contrast to percentages)
  for( i in length(out425[,1]):1 ) { #row
    out425$Pb[i]    = tb.m$black[i]  / tb.m$total[i];  # full blacks
    out425$R0.5[i]  = tb.m$B0.5[i]   / ( tb.m$total[i] - 2 * tb.m$black[i] );  # half blacks
  }

###generate plots
 tokens = unlist( strsplit(file, "/") );
 shortfile = tokens[2];
 #postscript( paste( "_", shortfile, "loh.042506.ps", sep=".") );

 par(mar=c(4,4,4,4));
 #plot the ratio
 y.max = max( out425[,2], na.rm=T ) * 1.2;
 y.min = min( out425[,2], na.rm=T ) ;
 cols = c("blue","black","green");
 if ( y.max == Inf ) { y.max = 7; }; plot(out425[,2] ~ out425$t,type='l', col=cols[1],
   main= paste( shortfile, " ", "Half vs Full Blacks" ),
   ylim=c(y.min, y.max), xlab="t (hours)",ylab="Half/Full Ratio");

  #overlay viability plot. 
  par(new=T);
  plot( out425$s ~ out425$t,type='l',xlab="",ylab="",axes=F,lty=2); 	

 #plot the percentage and risks
 y.max = max( out425[,3:4], na.rm=T ) * 1.2;
 y.min = min( out425[,3:4], na.rm=T ) ;
 par(new=T);
 plot(out425$Pb ~ out425$t, type='l',  col=cols[2],
   ylim=c(y.min, y.max), xlab="",ylab="",axes=F);
 lines( out425[,4] ~ out425$t, col=cols[3] );
 axis(4, at=pretty(c(y.min, y.max)) );
 mtext("Percentage", side=4, line=3);

 #legend
 labels = c("Half/Full", "P(full)", "R(1/2)", "viability");
 ltypes = c(1,1,1,2);
 legend( 0, y.max * 0.8, labels, col=c(cols,"black"), lty = ltypes);

 #dev.off();
} ###file loop####

dev.off();
# q("no");

