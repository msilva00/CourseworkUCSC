# Diana Gerardo
# Homework 5

library(MASS)
library(mvtnorm)
options(scipen = 999)
set.seed(859)

#########################################################################################
######################################### Part a ########################################
#########################################################################################
pois <- glm(faults ~ length, family=poisson(link='log'), data=fabric, x=T)

y <- fabric$faults
x <- fabric$length 

prior <- function(para){ return(log(1)) }

post <- function(para)
{
  return(exp(prior(para) + sum(y*(para[1]+para[2]*x) - exp(para[1]+para[2]*x))))
  
}

param <- matrix(NA, nrow=2, ncol=200000)
param[,1] <- c(0,0)
#param[,1] <- c(0.9717506,0.0019297)

new.param <- matrix(NA, nrow=2, ncol=200000)
new.param[,1] <- c(0,0)

J <- vcov(pois)
d <- 1


for(i in 2:200000)
{
  new.param[,i] <- rmvnorm(1, param[,i-1], d*J)
  
  q <- min(1, post(new.param[,i])/post(param[,i-1]))
  
  if(runif(1) < q)
  {param[,i] <- new.param[,i]}
  
  else
  {param[,i] <- param[,i-1]}
   
}


#Remove first 100000 iters
T <- 190000
param <- param[,-(1:T)] 


# Trace Plots and mean estimate
apply(param,1,quantile, c(0.025,0.5,0.975))

plot.ts(param[1,], ylab=expression(beta[1]),xlab="index")
mean(param[1,])
#0.958343

plot.ts(param[2,], ylab=expression(beta[2]), xlab="index")
mean(param[2,])
#0.001948863


# Density Plots
plot(density(param[1,]),
     xlab=expression(beta[1]),
     main=expression(beta[1]~" Distribution"))


plot(density(param[2,]),
     xlab=expression(beta[2]),
     main=expression(beta[2]~" Distribution"))


##################################
## estimated and interval plots ##
##################################

x.grid <- seq(min(x), max(x), length=100)

pi <- function(x,para){
  return(exp(para[1,]+para[2,]*x))
}

quantile <- matrix(NA,100,3)
for(i in 1:100){
  new.pi <- pi(x.grid[i],param)
  quantile[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
}

plot(x.grid,quantile[,1], type = "l", xlab="length of roll", 
     ylab="response mean", col="white",
     main="point and 95% interval estimates", ylim=c(0,20))


polygon(c(rev(x.grid), x.grid), 
        c(rev(quantile[ ,3]), quantile[ ,1]), 
        col =rgb(0,1,0,alpha=0.3), border = NA)

lines(x.grid,quantile[,2], col="darkgreen")

points(fabric$faults ~ fabric$length, pch=1, col="darkgreen")


# Posterior predictive distribution
postpred <- apply(param,2,function(theta){rpois(32,exp(theta[1]+theta[2]*x))})

plot(density(y),lwd=3, col="white", main="Orig. Data vs Post. Pred. y")
for(i in 1:10)
{
  lines(density(postpred[,i]), col="darkgreen")
}
lines(density(y),lwd=3)


# avg of all the post preds
for(i in 1:32)
{
  newpostpred[i] <- mean(postpred[i,])
}

# Orig Data vs Avg Post pred
plot(density(y),lwd=3, col="white", main="Orig. Data vs Avg. Post. Pred.", 
     xlab="faults",ylim=c(0,.13))
lines(density(newpostpred), col="darkgreen",lwd=2)
lines(density(y),lwd=3)



pred.sum = apply(postpred,1,quantile,c(.025,.975))
par(mfrow=c(1,1))
ind=1:length(y)
matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col="darkgreen", xlab="INDEX",ylab="faults",
        main="95% Posterior Predictive Interval Band ", ylim=c(0,29))
points(y, pch=19) 
identify(y ,n =3, labels=y)

### Residuals
res <- matrix(NA,ncol=32,nrow=2)

for(i in 1:32)
{
  res[1,i] <- quantile(y[i] - postpred[i,],.025)
  res[2,i] <- quantile(y[i]-postpred[i,],.975)
}

matplot(rbind(ind,ind), res, type="l", lty=1, col="darkgreen",
        main="Dist of Posterior Predictive Residuals", xlab="index",
        ylab="residuals")
abline(h=0)

#########################################################################################
######################################### Part b ########################################
#########################################################################################

set.seed(859)
pois <- glm(faults ~ length, family=poisson(link='log'), data=fabric, x=T)
y <- fabric$faults
x <- fabric$length 

# Priors
prior <- function(para){ return(log(1)) }
lam.prior <- function(lam){return(log(1/(1+lam)^2))}


# gamma : exp(beta1 + beta2*x)
gam <- function(betas)
{return(exp(betas[1]+betas[2]*x))}


# initializing
betas <- matrix(NA, nrow=2, ncol=50000)
betas[,1] <- c(0.9717506,0.0019297)


lambdas <- c()
lambdas[1] <- 1


mu_i <- matrix(NA, nrow=32, ncol=50000)
mu_i[,1] <- rnorm(32,0,1000)

J <- vcov(pois)
d <- 1

#power for prior analysis
p = 2


# k functions
# you get them by calulating the full conditionals of gamma_i and lambda

b.k <- function(lam,beta, mu)
{ return(sum(-lam*(beta[1]+beta[2]*x)-(mu*lam)/exp(beta[1]+beta[2]*x))) }

# this k function needs to have a +log(lambda) up in there outside the sum
# and i don't know where it comes from
# people keep saying from some Jacobian but idk what @__@
lam.k <- function(lam, beta, mu)
{
  return(-p*log(1+lam)+log(lam)+sum(lam*log(mu)-(mu*lam)/exp(beta[1]+beta[2]*x)+lam*log(lam)
                                    -lgamma(lam)-lam*(beta[1]+beta[2]*x)))
}



### Gibbin / MCMC ###
for(i in 2:50000)
{
  mu_i[,i] <- rgamma(32, y + lambdas[i-1], lambdas[i-1]/gam(betas[,i-1]) + 1)
  
  # M-H for betas
  new.betas <- rmvnorm(1, betas[,i-1], d*J)
  q1 <- min(1,exp(b.k(lambdas[i-1],new.betas,mu_i[,i])-b.k(lambdas[i-1], betas[,i-1], mu_i[,i])))
  if(runif(1) < q1)
  {betas[,i] <- new.betas}
  
  else
  {betas[,i] <- betas[,i-1]}
  
  
  #M-H for lambdas
  new.lambdas <- exp(rnorm(1, log(lambdas[i-1]),log(100)))
  q2 <- min(1,exp(lam.k(new.lambdas,betas[,i-1],mu_i[,i])-lam.k(lambdas[i-1], betas[,i-1], mu_i[,i])))
  if(runif(1) < q2)
  {lambdas[i] <- new.lambdas}
  
  else
  {lambdas[i] <- lambdas[i-1]}

}

# Remove first 10000 iters
T <- 10000
betas <- betas[,-(1:T)] 
lambdas <- lambdas[-(1:T)]
mu_i <- mu_i[,-(1:T)]

# Iteration and Density plots
plot.ts(betas[1,], xlab="iteration", ylab=expression(beta[1]))
plot(density(betas[1,]), xlab=expression(beta[1]), main="Distribution of "~beta[1])
mean(betas[1,])

plot.ts(betas[2,], xlab="iteration", ylab=expression(beta[2]))
plot(density(betas[2,]),xlab=expression(beta[2]), main="Distribution of "~beta[2])
mean(betas[2,])

# as "p" increases, lambda begins to follow a normal distribution
# it also affects the posterior predicive distribution 
plot.ts(lambdas, xlab="iteration", ylab=expression(lambda))
plot(density(lambdas),xlab=expression(lambda), main="Distribution of "~lambda)
mean(lambdas)
median(lambdas)


##########################################
## estimated and interval plots for (b) ##
##########################################
x.grid <- seq(min(x), max(x), length=100)

pi2 <- function(x,para){
  return(exp(para[1,]+para[2,]*x))
}

quantile2 <- matrix(NA,100,3)
for(i in 1:100){
  new.pi2 <- pi2(x.grid[i],betas)
  quantile2[i,] <- quantile(new.pi2,c(0.025,0.5,0.975))
}

plot(x.grid,quantile2[,1], type = "l", xlab="length of roll", 
     ylab="response mean", col="white",
     main="point and 95% interval estimates", ylim=c(0,20))


polygon(c(rev(x.grid), x.grid), 
        c(rev(quantile2[ ,3]), quantile2[ ,1]), 
        col =rgb(0,0,1,alpha=0.3), border = NA)

lines(x.grid,quantile2[,2], col="blue")

points(fabric$faults ~ fabric$length, pch=1, col="blue")
# note that the interval band is larger compared to part (a)



#######################################################################################
#################### Posterior Predictive Distribution with observed x ################
#######################################################################################

postpred2 <- apply(mu_i,2,function(theta){rpois(32,theta)})

# Orig. Data vs first 30 post pred dist
plot(density(y),lwd=3, col="white", main="Orig. Data vs Post. Pred. y")
for(i in 1:30)
{
  lines(density(postpred2[,i]), col="blue")
}
lines(density(y),lwd=3)

# Post Pred Check 95 interval
pred.sum2 = apply(postpred2,1,quantile,c(.025,.975))
par(mfrow=c(1,1))
ind=1:length(y)
matplot(rbind(ind,ind),pred.sum2,type="l",lty=1,col="blue", xlab="INDEX",ylab="faults",
        main="95% Posterior Predictive Interval Band ", ylim=c(0,40))
points(y, pch=19)  
identify(y ,n =2, labels=ind)

### Residuals
res2 <- matrix(NA,ncol=32,nrow=2)

for(i in 1:32)
{
  res2[1,i] <- quantile(y[i] - postpred2[i,], .025)
  res2[2,i] <- quantile(y[i]-postpred2[i,],.975)
}

matplot(rbind(ind,ind), res2, type="l", lty=1, col="blue",
        main="Dist of Posterior Predictive Residuals", xlab="index",
        ylab="residuals")
abline(h=0)



########################################################################################
######################## Posterior Predictive Distribution with NEW x ##################
########################################################################################

# you need a new x
# where x is fabric$length
x0 <- sample(min(x):max(x),32, replace = T)

# i think you also need new y's ???
# where y here stand for new fault values
y0 <- sample(min(y):max(y),32,replace = T)

# new gam fucntion for your new x
gam0 <- function(betas)
{ return(exp(betas[1]+betas[2]*x0)) }


# you need new mus
# lambdas and betas are your posterior samples
new.mu <- matrix(NA, nrow=32, ncol=40000)
for(i in 1:40000)
{
  new.mu[,i] <- rgamma(32, y + lambdas[i], lambdas[i]/gam(betas[,i]) + 1)  
}

# New posterior predictive
postpred2.2 <- apply(new.mu,2,function(theta){rpois(32,theta)})

### Residuals 
res2.2 <- matrix(NA,ncol=32,nrow=2)

for(i in 1:32)
{
  res2.2[1,i] <- quantile(y[i] - postpred2.2[i,], .025)
  res2.2[2,i] <- quantile(y[i]-postpred2.2[i,],.975)
}

matplot(rbind(ind,ind), res2.2, type="l", lty=1, col="blue",
        main="Dist of Posterior Predictive Residuals", xlab="index",
        ylab="residuals")
abline(h=0)


# Orig. Data vs first 30 post pred dist
plot(density(y),lwd=3, col="white", main="Orig. Data vs Post. Pred. y")
for(i in 1:10)
{
  lines(density(postpred2.2[,i]), col="blue")
}
lines(density(y),lwd=3)


# Post Pred Check 95 interval
pred.sum2.2 = apply(postpred2.2,1,quantile,c(.025,.975))
par(mfrow=c(1,1))
ind=1:length(y)
matplot(rbind(ind,ind),pred.sum2.2,type="l",lty=1,col="blue", xlab="INDEX",ylab="faults",
        main="95% Posterior Predictive Interval Band ", ylim=c(0,40))
points(y, pch=19)  
#identify(y ,n =2, labels=ind)


# avg of all the post preds
for(i in 1:32)
{
  newpostpred2[i] <- mean(postpred2.2[i,])
}

# Orig Data vs Avg Post pred
plot(density(y),lwd=3, col="white", main="Orig. Data vs Avg Post. Pred. ", ylim=c(0,.13))
lines(density(newpostpred2), lwd=2, col="blue")
lines(density(y),lwd=3)




#########################################################################################
######################################### Part c ########################################
#########################################################################################

## G & G

## Bayesian Poisson GLM
print(gg.m1 <- sum((y - apply(postpred,1,mean))^2) + sum(apply(postpred, 1, var)))
## Bayesian Poisson Hierarchical GLM
print(gg.m2 <- sum((y - apply(postpred2,1,mean))^2) + sum(apply(postpred2, 1, var)))


## Deviance Information Criterion

lpoislike <- function(y, mu){
  val = sum(dpois(y, mu,log=T))
  return(val)
}

### DIC M1 ###


