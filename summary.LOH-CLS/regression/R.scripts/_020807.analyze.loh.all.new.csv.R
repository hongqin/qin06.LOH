
rm(list=ls())

loh = read.table( "loh.new.csv", header=T, sep="\t");

loh$L0 = loh$half.0.raw / loh$b0

attach(loh) 

summary(lm(loh$Tc ~ loh$TLmax))
#Residual standard error: 2.042 on 37 degrees of freedom
#Multiple R-Squared: 0.3383,     Adjusted R-squared: 0.3204
#F-statistic: 18.92 on 1 and 37 DF,  p-value: 0.0001030

summary(lm(loh$Tdmax ~ loh$TLmax))
#Multiple R-Squared: 0.3726,     Adjusted R-squared: 0.3557
#F-statistic: 21.98 on 1 and 37 DF,  p-value: 3.684e-05

summary( lm( Lmax ~ TLmax + b0 + b.min + b.max + Tdmax + Tbmax + Tmmax + Tg + Tc))  #none

save.image("020807.loh.new.csv.RData");

