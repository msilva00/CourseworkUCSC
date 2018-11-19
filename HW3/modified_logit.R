################################################
## Monte Carlo Markov Chain implemented using the 
##  Metropolis Hastings Algorithm

################################################

############# Setting Up Workspace #############
rm(list = ls())
options(scipen = 999) ## Removing scientific notation
getwd()
## set.seed(1234) ## For testing purposes



## This produced a data frame with a binary dependent variable in the first column 
##  and the dosage level in the second
# data.frame(dead, dosage)
beetles = read.csv("beetles2.csv", header = T, sep = ',')
dead = beetles$dead
dosage = beetles$dosages

############# Running the Markov Chain ############

## Number of Steps for MCMC
N = 500000

## Two vectors storing the beta coefficients
b0 <- numeric()
b1 <- numeric()
alpha <- numeric()

## Set starting values for b0 and b1
b0[1] = 0
b1[1] = 0
alpha[1] = 1


###################################################################
#### This is the part we modify (i.e. logit to modified logit) ####
loglikelihood <- function(b0, b1, alpha){
  ## Calculate p through link function
  p = exp(alpha*(b0+b1*dosage))/((1 + exp(b0 + b1*dosage))^alpha)
  # print(p)
  ## sum dbinom at each level of dosage
  sum(dbinom(dead,1, p, log=TRUE))
}
###################################################################

## Generating N draws from random uniform for acceptance/rejection
## This is more computationally efficient than drawing from runif in every realization
uniforms <- runif(N, 0, 1)

## Running the algorithm for N realizations
for (i in 2:N) {
  ## Draw from the candidate, which is a normal dist. centered around previously accepted value
  ## SDs were approximated using trial and error (they were tuned)
  b0.prime <- rnorm(1, b0[i-1], 1)  
  b1.prime <- rnorm(1, b1[i-1], 1)
  alpha.prime <- rgamma(1, alpha[i-1], 2)
  
  ## Setting prior densities
  ## I didn't have prior knowledge of beetle deaths, so I used a flat prior
  ## i.e. a normal distribution with a very large standard deviation
  current.prior.dens = dnorm(b0.prime,0,100000,log = TRUE) + dnorm(b1.prime,0,100000,log=TRUE) + dgamma(alpha.prime, 1,1)
  prime.prior.dens = dnorm(b0[i-1],0,100000, log = TRUE) + dnorm(b1[i-1],0,100000, log = TRUE) + dgamma(alpha[i-1], 1,1)
  
  ## Calculating log likelihood for both the proposed and current values of b0 and b1
  current.loglik = loglikelihood(b0[i-1], b1[i-1], alpha[i-1])
  prime.loglik = loglikelihood(b0.prime, b1.prime, alpha.prime)
  
  ## Summing the prior and likelihood
  current.target.dens = current.loglik + current.prior.dens
  prime.target.dens = prime.loglik + prime.prior.dens
  
  ## Alpha denotes the acceptance probability 
  acceptance =  exp(prime.target.dens - current.target.dens)
  
  ## Comparing alpha against a draw from U[0,1]
  if (uniforms[i] < acceptance) {
    ## Updating the chain if proposed value is accepted
    b0[i] = b0.prime
    b1[i] = b1.prime
    alpha[i] = alpha.prime
  } else {
    ## Updating the chain if chain proposed value is not accepted (chain doesn't move)
    b0[i] = b0[i-1]
    b1[i] = b1[i-1]
    alpha[i] = alpha[i-1]
  }
}

## Create a data frame which holds the chains for both parameters
results.vec <- data.frame(b0, b1, alpha)

layout(matrix(1:3, ncol = 1), widths = 1,
       heights = c(5,5,5), respect = FALSE)
par(mar=c(2.5, 4, 1, 1),oma=c(0,1,2,0))
mtext("Metropolis-Hastings Traceplot")
plot(results.vec[, 1], type = "l", ylab = expression(beta[1]))
plot(results.vec[, 2], type = "l", ylab = expression(beta[2]))
plot(results.vec[, 3], type = "l",ylab = expression(alpha))


plot(1, type = 'n', axes = FALSE, bty = 'n', ylab = '')
legend('left', , c("X","Y"), bty="n", horiz=T, cex=1.5, col=c("red1","darkblue"), text.col=c("red1","darkblue"), pch=c(1,3), lty=c(2,3), x.intersp=0.4,adj=0.2)
par(mar=c(0, 4, 2, 1), bty = 'o')
plot(a1, type="b", ylim=c(0,14.5), xlab="Time (secs)", ylab="", xaxt = 'n', cex.axis=1.4, cex.lab=1.3,cex=1.2, lwd=2.5, col="red1", lty=2, pch=1, main="A")
lines(a2,type="b",pch=3,lty=3,col="darkblue",lwd=2.5,cex=1.2)
par(xpd=T)
plot(b1, type="b", ylim=c(0,14.5), xlab="Time (secs)", ylab="", xaxt = 'n', cex.axis=1.4, cex.lab=1.3,cex=1.2,lwd=2.5,col="red1",lty=2,pch=1, main="B")
lines(b2,type="b",pch=3,lty=3,col="darkblue",lwd=2.5,cex=1.2)
plot(c1, type="b", ylim=c(0,14.5), xlab="Time (secs)", ylab="", xaxt = 'n', cex.axis=1.4, cex.lab=1.3,cex=1.2,lwd=2.5,col="red1",lty=2,pch=1, main="C")
lines(c2,type="b",pch=3,lty=3,col="darkblue",lwd=2.5,cex=1.2)
par(mar=c(4, 4, 2, 1))
plot(d1/1000, type="b", ylim=c(0,14.5), xlab="Time (secs)", ylab="", xaxt = 'n', cex.axis=1.4, cex.lab=1.3,cex=1.2,lwd=2.5,col="red1",lty=2,pch=1, main="D")
lines(d2,type="b",pch=3,lty=3,col="darkblue",lwd=2.5,cex=1.2)
mtext("Price", side=2, at=40,line=2.5,cex=1.1)
axis(1, 1:10, cex.axis = 1.4)

## Plotting each of the parameters
par(mfrow=c(2,1))


## From the plot, the approximate burn-in period for the intercept is 35
## Burn-in for b1 is approximately 70

## These tend to change from sample to sample because of randomization
b0.burnin = 100000
b1.burnin = 100000
alpha.burnin = 100000

plot(results.vec[-(1:b0.burnin), 1], type = "l")
plot(results.vec[-(1:b1.burnin), 2], type = "l")
plot(results.vec[-(1:alpha.burnin), 3], type = "l")

## Discarding burn-in samples
b0.new <- b0[-(1:b0.burnin)]
b1.new <- b1[-(1:b1.burnin)]
alpha.new<- alpha[-(1:alpha.burnin)]
# ## Checking autocorrelation using ACF
# acf(b0.new, lag.max = 500)
# acf(b1.new, lag.max = 500)

## From ACF on both parameters, to thin both parameters, I will take roughly every 110th sample for the b0
## and every 120th sample for b1. This is because the autocorrelation dies out after those many lags.

## These tend to change from sample to sample because of randomization
b0.thin = 110 
b1.thin = 120
alpha.thin = 100

## Thinning the chains
b0.new <- b0.new[seq(1, length(b0.new), b0.thin)]
b1.new <- b1.new[seq(1, length(b1.new), b1.thin)]
alpha.new <- alpha.new[seq(1, length(alpha.new), alpha.thin)]

## Displaying summary stats
summary(b0.new)
summary(b1.new)
summary(alpha.new)

par(mfr)
## Displaying thinned histograms for each parameter
hist(b0.new, freq = FALSE, breaks = 100)
hist(b1.new, freq = FALSE, breaks = 100)
hist(alpha.new, freq = FALSE, breaks = 100)
