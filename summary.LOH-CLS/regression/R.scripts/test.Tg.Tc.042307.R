N = 110;

Tc = rnorm( N, mean=4, sd=1);

Tg = 1.2 * Tc + 0.8 + rnorm(N, sd=1.7);
summary( lm( Tg ~ Tc ) )

Tg.vs.Tc = Tg/ Tc

arls =  12.6 * Tg.vs.Tc + 14.2 + rnorm(N, sd=3.1 )

plot( arls ~ Tg )

plot( arls ~ Tc )

summary( lm(arls ~ Tg) );

summary( lm(arls ~ Tc) );

# small samle size are more influenced by noise.
# if RM11 is removed, there is a negative correlation between 
# Arls and ClS. 

Tg = 1.2 * Tc + 0.8
ARLS = 12.6 Tg.vs.Tc + 14.2


###########3
A = 1.5
B = 5
x = rnorm(100)
y = A * x + B + rnorm(100)
summary( lm( y ~ x) )





