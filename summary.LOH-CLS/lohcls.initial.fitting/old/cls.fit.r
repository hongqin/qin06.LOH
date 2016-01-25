library(nlme)
library(survival)

cls = read.table ("cls.curve.fit.tab", header=T);
cls[ , c(1,3,5,7)] = cls[ , c(1,3,5,7)] / 24

logistical.viability <- function( v, w, t ) { ret <- 1 /( 1 + ( t / v )^ w );  }

cls[ ,2] = cls[ ,2] / 1.1
cls[1,2]=1

xlim = c(0, 20);
ylim = c(0,1.1);




######## fitting RM11
 x = seq(0,20, 0.02);
 s <- cls[ ,4]; s = s[! is.na(s)];
 t <-  cls[ ,3]; t =  t[! is.na(t)];
 fm.log <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 10, w = 4 )  );
 summary( fm.log  );

v.fit <- fm.log$coefficients[1]; v.fit
w.fit <- fm.log$coefficients[2]; w.fit

s.fit <- logistical.viability( v.fit, w.fit, x ) ;

plot ( cls[2:17,4]~cls[2:17,3], type='p' , xlim=xlim, ylim=ylim, ylab="viability", xlab="days", pch=24, cex=0.75);
#points( cls[,4]~cls[,3], pch=24, cex=0.75);
lines( s.fit ~ x );


######## fitting SGU
 x = seq(0,15, 0.005);

 s <- cls[ ,2]; s = s[! is.na(s)];
 t <-  cls[ ,1]; t =  t[! is.na(t)];
 fm.log <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 10, w = 4 )  );
 summary( fm.log  );

v.fit <- fm.log$coefficients[1]; v.fit
w.fit <- fm.log$coefficients[2]; w.fit

s.fit <- logistical.viability( v.fit, w.fit, x ) ;

#plot ( cls[,2]~cls[,1], type='p' , xlim=xlim, ylim=ylim, ylab="viability", xlab="days");
points( cls[2:17,2]~cls[2:17,1], pch=25, cex=0.75);
lines( s.fit ~ x );

######## fitting YPS128
 x = seq(0,11, 0.02);

 #cls[ , 6] 
 s <- cls[ ,6]; s = s[! is.na(s)];
 t <-  cls[ ,5]; t =  t[! is.na(t)];
 fm.log <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 10, w = 4 )  );
 summary( fm.log  );

v.fit <- fm.log$coefficients[1]; v.fit
w.fit <- fm.log$coefficients[2]; w.fit

s.fit <- logistical.viability( v.fit, w.fit, x ) ;

#plot ( cls[,6]~cls[,5], type='p' , xlim=xlim, ylim=ylim, ylab="viability", xlab="days");
points( cls[2:17,6]~cls[2:17,5], pch=21);
lines( s.fit ~ x );


######create legends
shape = c( 24,21,25,);
strains = c("RM11", "YPS128","SGU57", );
legend( 15, 1.0, strains, pch=shape, );
