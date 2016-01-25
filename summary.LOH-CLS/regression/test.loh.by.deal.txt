 library(deal);

 sub = tb[,c("ARLS","b.max","Tg.vs.Tc")]
 fit = network(sub);

 fit = network(sub, specify=T);  #disable ARLS~b.max, but this is also problematic. 

 fit.prior <- jointprior(fit, 11)
 
 fit       <- getnetwork(learn(fit, sub, fit.prior))

 plot(fit); # this gives only 3 vertices

 hisc      <- autosearch(fit, sub, fit.prior, trace=T)

 plot(getnetwork(hisc))

 hisc2      <- autosearch(fit, sub, fit.prior,trace=FALSE,removecycles=TRUE) # slower

 plot( getnetwork(hisc2) );

Read through DEAL's paper and manual. It seems that DEAL does not use d-separation in the topology search. 


#############################


 summary( lm( rats$W1 ~ rats$W2 ) );
#R2=0.51, p=8.8E-5

 summary( lm( rats$W1 ~ rats$Drug ) );
#R2=0.73, p=9.4E-7

 summary( lm( rats$W2 ~ rats$Drug ) );
#R2=0.298, p=0.0983


 summary( lm( rats$W2 ~ rats$W1 + rats$Drug ) );
#R2=0.68, p=3.34E-5
#for W1, p=2.16E-5
#p=0.36 for Drug2
#p=0.017 for Drug1
#Hmm, this seems to suggest that Drug and W2 is d-separated by W2.
#but, how is the direction determined?


 m1 = lm( rats$W2 ~ rats$Drug);
 m2 = lm( rats$W1 ~ rats$Drug);
 summary( lm( m1$residuals ~ m2$residuals) );
#R2=0.60, p=8.2E-6


 m1 = lm( rats$W2 ~ rats$Drug);
 m2 = lm( rats$W1 ~ rats$Drug);
 summary( lm( m1$residuals ~ m2$residuals) );


###################
data(rats)
fit       <- network(rats)
fit.prior <- jointprior(fit,12)
fit       <- getnetwork(learn(fit,rats,fit.prior))
plot(fit); # this gives only 4 vertices
hisc      <- autosearch(fit,rats,fit.prior,trace=FALSE)
plot(getnetwork(hisc))

