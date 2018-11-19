
# Data 2: Revisit the Beetle Mortality Data (Dobson)
# y: Number of beetles killed after fives hours exposure to
#    gaseous carbon disulphide
# n: Number of beetles 
# x: Dose

beetle<-as.data.frame(matrix(c(1.6907, 59, 6,
                               1.7242, 60, 13,
                               1.7552, 62, 18,
                               1.7842, 56, 28,
                               1.8113, 63, 52,
                               1.8369, 59, 53,
                               1.861, 62, 61,
                               1.8839, 60, 60), ncol = 3, byrow=T))
colnames(beetle)<-c("dose", "n", "kill")

#----------------------------------------
# 1. Residuals and residual plots
#----------------------------------------


# 1.2 Beetle data (counts)

plot(beetle$dose, beetle$kill/beetle$n, pch=19, cex=1.5, main="Beetle Mortality", ylab="Mortality Rate")
##### Logit link ######
beetle.logit<-glm(kill/n~dose, family=binomial(link=logit), weight=n,
                  data=beetle)

# Predictions and residuals
pi <- beetle.logit$fitted  
lp <- predict(beetle.logit, type="link")
yhat<-beetle$n*pi
e <- beetle$kill-yhat  # Oridinary residual
rp<- e/sqrt(beetle$n*pi*(1-pi))  # Pearson Residual
dev <- residuals(beetle.logit)  # Deviance residual
beetle.res1 <- cbind(beetle$kill,beetle$yhat1, pi,lp, e,rp,dev)
round(beetle.res1, 4)

# Prepare graphs
par(mfrow=c(2,2), las=1)
indi <- sort(pi,index.return=T)$ix
loe<-loess(rp ~ seq(1,8), degree=1)
plot(rp, type="p", pch=16, xlab="Index.", ylab="Pearson Residual", main="Logit Link")
lines(seq(1,8), loe$fitted[indi], type="l",col="red")
abline(h=0, col = "gray")

pi = seq(1,8)
loe<-loess(dev ~ pi, degree=1)
plot(pi, dev, type="p", pch=16, xlab="Index", ylab="Deviance Residual", main="Logit Link")
lines(pi[indi], loe$fitted[indi], type="l",col="red")
abline(h=0, col = "gray")

##### probit link ######
beetle.probit<-glm(kill/n~dose, family=binomial(link=probit), weight=n, data=beetle)
pi <- beetle.probit$fitted  # Estimated probability
yhat<-beetle$n*pi
e <- beetle$kill-yhat  # Oridinary residual
rp<- e/sqrt(beetle$n*pi*(1-pi))  # Pearson Residual
dev <- residuals(beetle.probit)  # Deviance residual
indi <- sort(pi,index.return=T)$ix
pi = seq(1,8)
loe<-loess(rp ~ pi, degree=1)
plot(pi, rp, type="p", pch=16, xlab="Index", ylab="Pearson Residual", main="Probit Link")
lines(pi[indi], loe$fitted[indi], type="l", col = "red")
abline(h=0, col = "gray")

loe<-loess(dev ~ pi, degree=1)
plot(pi, dev, type="p", pch=16, xlab="Index", ylab="Deviance Residual", main="Probit Link")
lines(pi[indi], loe$fitted[indi], type="l", col = "red")
abline(h=0, col = "gray")

##### c-log-log link ######
beetle.clog<-glm(kill/n~dose, family=binomial(link=cloglog), weight=n, data=beetle)
pi <- beetle.clog$fitted  # Estimated probability
yhat<-beetle$n*pi
e <- beetle$kill-yhat  # Oridinary residual
rp<- e/sqrt(beetle$n*pi*(1-pi))  # Pearson Residual
dev <- residuals(beetle.clog)  # Deviance residual
indi <- sort(pi,index.return=T)$ix

pi = seq(1,8)
loe<-loess(rp ~ pi, degree=1)
plot(pi, rp, type="p", pch=16, xlab="Index", ylab="Pearson Residual", main="C-LogLog Link")
lines(pi[indi], loe$fitted[indi], type="l", col = "red")
abline(h=0,col="gray")

loe<-loess(dev ~ pi, degree=1)
plot(pi, dev, type="p", pch=16, xlab="Index", ylab="Deviance Residual", main="C-LogLog Link")
lines(pi[indi], loe$fitted[indi], type="l", col = "red")
abline(h=0, col = "gray")

#----------------------------------------
# 3. Goodness of fit tests
#----------------------------------------

# 3.1. Hosmer-Lemeshow test (Appropriate for 0-1 responses)
# The following function is addopted from
# http://www.math.mcmaster.ca/peter/s4f03/s4f03_0607/rochl.html
hosmerlem<-function(y, yhat, g = 10)
{
  cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0, 1, 1/g)), include.lowest = T)
  obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2/expect)
  P <- 1 - pchisq(chisq, g - 2)
  c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
} 
# 3.2. Pearson Chi-square (counts data only)
pgof<-function(n, y, pihat, p)
{
  Oy<-y
  On<-n-y
  Ey<-n*pihat
  En<-n*(1-pihat)
  chisq<-sum((Oy-Ey)^2/Ey) + sum((On-En)^2/En)
  pvalue <- 1 - pchisq(chisq, length(n)-p)
  c("X^2" = chisq, Df = length(n)-p, "P(>Chi)" = pvalue)
}
# Use Beetle Mortality data again
logit_P2 = pgof(beetle$n, beetle$kill, beetle.logit$fitted, 2)
probit_P2 = pgof(beetle$n, beetle$kill, beetle.probit$fitted, 2)
cloglog_P2 = pgof(beetle$n, beetle$kill, beetle.clog$fitted, 2)
Chi_sq=as.data.frame(rbind(logit_P2,probit_P2,cloglog_P2))

# 3.3. Deviance test (counts data only)
c(beetle.logit$dev, sum(residuals(beetle.logit)^2))  # just to show they are the same
1-pchisq(beetle.logit$dev, df=length(beetle$dose)-2) # why -2?
1-pchisq(beetle.clog$dev, df=length(beetle$dose)-2) # why -2?

#------------------------------------------------
# This is the end of script 4
#------------------------------------------------