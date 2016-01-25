# output basal loh and max loh

rm(list=ls());

file = "INPUTFILE";
#file = "normalized.data/norm.M28.022006.txt";

tb.m = read.table( file, header=T, sep="\t");

# todo: remove NA and zeros from tb.m
# output to out3,
# columns in out3:
 out3 = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= 7 ) );
 names( out3 ) = c( "t", "half.over.black", "quarter.over.half", "Pb", "R0.5", "R0.25", "R2" );

 out3[,1] = tb.m[,1]; # "t"

 #calculate the ratios
  for( i in length(out3[,1]):1 ) { #row
    out3[i,2] = tb.m$B0.5[i]  / tb.m$black[i];   
    out3[i,3] = tb.m$B0.25[i] / tb.m$B0.5 [i]; 
  }

 #calculate the risks (in contrast to percentages)
  for( i in length(out3[,1]):1 ) { #row
    out3$Pb[i]    = tb.m$black[i]  / tb.m$total[i];  # full blacks
    out3$R0.5[i]  = tb.m$B0.5[i]   / ( tb.m$total[i] - 2 * tb.m$black[i] );  # half blacks
    out3$R0.25[i] = tb.m$B0.25[i]  / ( tb.m$total[i] - 2 * tb.m$black[i] );  # quarter blacks
    out3$R2[i]    = tb.m$B2[i]     / ( tb.m$total[i] - 2 * tb.m$black[i] );  # quarter blacks
  }

  tokens = unlist( strsplit(file, "/") );
  shortfile = tokens[2];
  outfile = paste("loh.risk/loh.", shortfile, sep="");
  write.table( out3, outfile, quote=F, row.names=F, sep="\t");

###generate plots
 postscript( paste( "_", shortfile, "loh.risk.ps", sep=".") );

 #plot the ratio
 y.max = max( out3[, 2:3], na.rm=T ) * 1.2;
 y.min = min( out3[, 2:3], na.rm=T ) ;
 cols = c("black","green");
 if ( y.max == Inf ) { y.max = 20 };

 plot(out3[,2] ~ out3$t,type='l', col=cols[1],
   main= paste( shortfile, " ", "ratios" ),
   ylim=c(y.min, y.max), xlab="t",ylab="ratios");
 lines( out3[,3] ~ out3$t, col=cols[2]);
 labels = names( out3 );
 legend( 0, y.max , labels[2:3], col=cols, lty=1);

 #plot the percentage and risks
 y.max = max( out3[, 4:7], na.rm=T ) * 1.2;
 y.min = min( out3[, 4:7], na.rm=T ) ;
 cols = c("black","green","red","purple","blue");
 # par( yaxp=c(1e-4,20,1));
 plot(out3[,4] ~ out3$t,type='l',  col=cols[1],
   main= paste( shortfile, " ", "risk" ),
   ylim=c(y.min, y.max), xlab="t",ylab="percentage and risk");
 count = 1;
 for ( i in 5:7 ) {
   count = count + 1;
   lines( out3[,i] ~ out3$t, col=cols[ count ] );
 }
 labels = names( out3 );
 legend( 0, y.max * 0.75, labels[4:7], col=cols, lty=1);
 dev.off();

# q("no");

