#logistic model of P(black)

rm( list=ls() );

#file = "INPUTFILE";
file = "normalized.data/norm.M28.022006.txt";

tb.m = read.table( file, header=T, sep="\t");

# todo: remove NA and zeros from tb.m
# output to out3,
# columns in out4.21:
 header = c("t","s","Pb","Rb","R0.5.vs.Rfull","R0.5","R0.25","R2", "Pb.fm");
 out4.21 = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out4.21 ) = header;

 out4.21[,1] = tb.m[,1]; # "t"

 #simple calculations
  for( i in length(out4.21[,1]):1 ) { #row
    out4.21$s[i]     = tb.m$total[i]  / tb.m$total[1]; #s
    out4.21$Pb[i]    = tb.m$black[i]  / tb.m$total[i];  # full blacks
    out4.21$R0.5[i]  = tb.m$B0.5[i]   / ( tb.m$total[i] - 2 * tb.m$black[i] );  # half blacks
    out4.21$R0.25[i] = tb.m$B0.25[i]  / ( tb.m$total[i] - 2 * tb.m$black[i] );  # quarter blacks
    out4.21$R2[i]    = tb.m$B2[i]     / ( tb.m$total[i] - 2 * tb.m$black[i] );  # quarter blacks
  }

 #fit percentage of blacks Pb (without the first time point?)
 library(nlme)
 attach(out4.21);
 a.start = max( Pb );
 b.start = a.start/Pb[1] - 1;
 c2 =  t[2] / log(b.start/(a.start/Pb[2]-1 ));
 c3 =  t[3] / log(b.start/(a.start/Pb[3]-1 ));
 c.start = 0.5 * c2 + 0.5 * c3;
# fm <- gnls( Pb ~ a / (1 + b * exp(-t/c) ), start= list(a=a.start, b=b.start,c=c.start));
 fm <- gnls( Pb ~ a / (1 + b /exp((t-tc)/c) ), start= list(a=a.start, b=b.start,c=c.start,tc=80));
 a.fm = fm$coefficients[1];
 b.fm = fm$coefficients[2];
 c.fm = fm$coefficients[3];
 tc.fm = fm$coefficients[4];
 Pb.fun = function( a, b, c, tc, t ) { ret <- a / (1 + b * exp(-(t-tc)/c) ) ; }
 Pb.fm = Pb.fun( a.fm, b.fm, c.fm, tc, t);

 a.start = max(Pb);
 c.start = Pb[1];
 b.start = 80;
 fm2 <- gnls( Pb ~ a * t/(t + b) + c, start= list(a=a.start, b=b.start,c=c.start));

