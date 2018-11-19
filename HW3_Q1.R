setwd("Desktop/")

#1a
# Data
beetles <- read.csv("beetles.txt", sep = '' , header = F)
names(beetles) = c("logDose", "numBeetles", "numDied")

beetles
dose = beetles$logDose
killed = beetles$numDied
number = beetles$numBeetles

beetles$prob <- beetles$numDied/beetles$numBeetles

# Exploratory Analysis

pairs(beetles[,c("logDose","prob")])

# Fit models
link <- c("logit","probit","cloglog")
models <- lapply(as.list(link), function(lnk){
  glm(prob~logDose, 
      weights=numBeetles, 
      family=binomial(link=lnk),
      data=beetles)
})

glm_logit =glm(prob~logDose, 
    weights=numBeetles, 
    family=binomial(link="logit"),
    data=beetle)
summary(glm_logit)
confint(glm_logit)
?confint
plot(glm_logit)

xgrid = seq(-3,4,100)

plot(fitted(glm_logit),(beetles$numDied/beetles$numBeetles),
     xlab="Fitted",ylab="Observed", main = "Logit Link")
abline(a=0,b=1)
?abline

sort(beetles$numDied/beetles$numBeetles)

glm_probit =glm(prob~logDose, 
               weights=numBeetles, 
               family=binomial(link="probit"),
               data=beetles)
summary(glm_probit)
confint(glm_probit)

plot(fitted(glm_probit),(beetles$numDied/beetles$numBeetles),
     xlab="Fitted",ylab="Observed", main = "Probit Link")
abline(a=0,b=1)


fit.cloglog =glm(prob~logDose, 
                weights=numBeetles, 
                family=binomial(link="cloglog"),
                data=beetles)
summary(fit.cloglog)
confint(fit.cloglog)

plot(fitted(fit.cloglog),(beetles$numDied/beetles$numBeetles),
     xlab="Fitted",ylab="Observed", main = "C-Log-Log Link")
abline(a=0,b=1)

predict=fitted(fit.cloglog) # predicted proportions
pearson.resid=resid(fit.cloglog,type="pearson")
pearson.resid.std=pearson.resid/sqrt(1-lm.influence(fit.cloglog)$hat)
matrix(cbind(dose,killed/number,predict,pearson.resid.std),ncol=4,
       dimnames=list(NULL,c("dose","proportion","fitted","std. pearson resid")))


names(models) <- sapply(models,function(m) m$family$link)

# Residual Analysis: mean squared deviance residuals
dev_resid <- lapply(models,function(x) resid(x,type="deviance"))
sapply(dev_resid,function(x) sqrt(mean(x^2)))

# Pred_Pi
pred_pi <- sapply(models,predict,type="response")
# pdf("../img/freqpreds.pdf")
plot(beetles$prob,pred_pi[,1],ylim=0:1,type='n',fg='grey',bty='n',
     ylab="pi hat",xlab="observed probability")
for (i in 1:3) 
  points(beetles$prob,pred_pi[,i],pch=i+1,cex=2,lwd=3,col='grey30')
 
abline(0,1,col='grey')
legend("bottomright",legend=link,pch=2:4,cex=2,
       col='grey', text.col='grey',box.col='grey')
dev.off()

## grid predictions
xx <- seq(1.65, 1.9, length=100)
pred_pi_grid<-sapply(models,predict,newdata=list(logDose=xx),type="response")

#1c
inv.mod.logit <- function(x, b, a) {
  eta <- c(x %*% b)
  (exp(eta) / (1+exp(eta))) ^ a
}
b.hat <- c(-113.625, 62.5)
a.hat <- .279
pred_mod_logit <- inv.mod.logit(cbind(1,xx), b.hat, a.hat)
# end of 1c

pdf("../img/freqcurves.pdf")
plot(xx,pred_pi_grid[,1],ylim=0:1,type='n',fg='grey',bty='n',
     ylab="pi hat grid",xlab="observed probability")
for (i in 1:3) 
  lines(xx,pred_pi_grid[,i],col=i+1,lwd=3)
abline(0,1,col='grey')
points(beetles$logDose,beetles$prob,pch=20,col='grey30',cex=2)
lines(xx,pred_mod_logit,col="pink",lwd=3)
legend("topleft",legend=c("Data",link,"modified logit"),cex=2,text.col=c("grey30",2:4,"pink"), bty='n')
dev.off()

p.hat <- inv.mod.logit(cbind(1,beetles$logDose), b.hat, a.hat)
m <- beetles$numBeetles
y.hat <- p.hat * m
y <- beetles$numDied
n <- length(y)
dev.new <- sign(y-y.hat) * sqrt(2*abs(y*log(y/y.hat) + (m-y)*log(1E-10+(m-y)/(m-y.hat))))


pdf("../img/resid2.pdf")
plot(dev.new,col="black",pch=1,cex=1,ylab="deviance residuals")
#plot(dev_resid$cloglog,pch=20,cex=3)
#points(dev.new,col="pink",pch=20,cex=3,ylab="deviance residuals")
abline(h=0)
dev.off()

#1d
logchoose <- function(m,y) lgamma(m+1) - lgamma(y+1) - lgamma(m-y+1)


K <- 3 # number of parameters
aic <- c(sapply(models,function(m) m$aic),-2 * sum(logchoose(m,y) + y*log(p.hat) + (m-y)*log(1-p.hat)) + 2*K)
names(aic)[4] <- "logit2"
bic <- aic + c(2,2,2,3) * (log(n)-2)

aic
bic

