ARLS = tb$ARLS

m1 = lm( ARLS ~ Tg.vs.Tc );
m2 = lm( ARLS ~ b.max );
m3 = lm( Tg.vs.Tc ~ b.max );

summary(m1)
summary(m2);
summary(m3);
#all correlated, therefore, this is not a collider

Ma = lm( ARLS ~ Tg.vs.Tc + b.max );
Ma = lm( ARLS ~ b.max + Tg.vs.Tc );
summary(Ma);
# ARLS Indepedent of b.max | Tg.vs.Tc

Mb = lm( Tg.vs.Tc ~ ARLS + b.max );
summary(Mb);
# Tg.vs.Tc indepdent of b.max | ARLS

Mc = lm( b.max ~ Tg.vs.Tc + ARLS );
summary(Mc);

#conclusion:
# 