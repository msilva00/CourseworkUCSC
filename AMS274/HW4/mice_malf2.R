mice <- read.csv("mice.csv")
mice$p1 <- mice$dead/mice$numSubj
mice$n2 <- mice$numSubj - mice$dead
mice$p2 <- mice$malf/mice$n2
mice$p3 = mice$normal/mice$numSubjects



mod2 <- glm(p2 ~ conc, 
            weights=n2,
            family=binomial(link=logit), data=mice)
summary(mod2)


pi2 = predict(mod2,newdata=list(conc=newconc),type="response")


newconc <- seq(0, 500, length.out = 100)
p1 <- predict(mod2,newdata=list(conc=newconc),type="response")
pp <- predict(mod2,newdata=list(conc=newconc),se.fit = T)
linkinv1 <- family(mod2)$linkinv ## inverse-link function

pred0 <- pp$fit
pred <- linkinv1(pp$fit)
alpha <- 0.95
sc <- abs(qnorm((1-alpha)/2)) 
alpha2 <- 0.5
sc2 <- abs(qnorm((1-alpha2)/2))
pframe = cbind(pred0,pred)
pframe <- transform(pframe,
                    lwr=linkinv1(pred0-sc*pp$se.fit),
                    upr=linkinv1(pred0+sc*pp$se.fit),
                    lwr2=linkinv1(pred0-sc2*pp$se.fit),
                    upr2=linkinv1(pred0+sc2*pp$se.fit))
with(pframe,
     {
       plot(newconc,pred,ylim=c(0,1), col = "white", lty = 3)
       arrows(as.numeric(newconc),lwr,as.numeric(newconc),upr,
              angle=90,code=3,length=0.1, col = "red")
     })

plot(mice$conc, seq(0,1,length=length(mice$conc)),type='n', fg='grey', 
     xlab='concentration (mg/kg per day)', ylab='probability',
     las=1, main=expression(hat(pi[j])(x)),cex.main=1,ylim=c(0,1))
lines(newconc, pi2, col='green', lwd=1)
lines(newconc, linkinv1(pred0+sc*pp$se.fit), col = 'green', lwd = 1, lty = 6)
lines(newconc, linkinv1(pred0-sc*pp$se.fit), col = 'green', lwd = 1, lty = 6)





