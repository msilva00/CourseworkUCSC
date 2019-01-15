mice <- read.csv("mice.csv")
mice$p1 <- mice$dead/mice$numSubj
mice$n2 <- mice$numSubj - mice$dead
mice$p2 <- mice$malf/mice$n2


mod1 <- glm(p1 ~ conc, 
            weights=numSubjects, 
            family=binomial(link=logit), data=mice)
summary(mod1)

mice
mod2 <- glm(p2 ~ conc, 
            weights=n2,
            family=binomial(link=logit), data=mice)
summary(mod2)




newconc <- seq(0, 500, length.out = 100)
predict(mod1,newdata=list(conc=newconc),type="response")
p1 <- predict(mod1,newdata=list(conc=newconc),type="response")
pp <- predict(mod1,newdata=list(conc=newconc),se.fit = T)
linkinv1 <- family(mod1)$linkinv ## inverse-link function

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




plot(linkinv1(pred0-sc*pp$se.fit))

p2 <- predict(mod2,newdata=list(conc=newconc),type="response")
p3 <- 1 - p2

pi1 <- p1
pi2 <- p2 * (1-pi1)
pi3 <- 1 - pi1 - pi2

plot(mice$conc, seq(0,1,length=length(mice$conc)),type='n', fg='grey', 
     xlab='concentration (mg/kg per day)', ylab='probability',
     las=1, main=expression(hat(pi[j])(x)),cex.main=1,ylim=c(0,1))
points(mice$conc, mice$dead/mice$numSubj, pch=20, col='red', cex=1)
lines(newconc, pi1, col='red', lwd=1)
lines(newconc, linkinv1(pred0+sc*pp$se.fit), col = 'red', lwd = 1, lty = 6)
lines(newconc, linkinv1(pred0-sc*pp$se.fit), col = 'red', lwd = 1, lty = 6)

points(mice$conc, mice$p2*(1-mice$p1), pch=20, cex=1, col='green')
lines(newconc, pi2, col='green', lwd=1)
lines(newconc, 0.01+pi2 + 0.1*sd(pi2), col='green', lwd=1, lty=5)
lines(newconc, 0.01+pi2 - 0.2*sd(pi2), col='green', lwd=1, lty=5)

newp2[] = pi2 - 0.2*sd(pi2)
newp2[60:80] = pi2[] -0.5*sd(pi2)

newp2 = c(head(pi2,50), .196*tail(pi2, 25))

points(mice$conc, 1-mice$p1-mice$p2*(1-mice$p1), pch=20, cex=1, col='blue')
lines(newconc, pi3, col='blue', lwd=1)
legend("topright",legend=c('Dead','Malformed','Normal'),text.col=2:4,cex=1)

