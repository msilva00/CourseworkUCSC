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

## Set starting values for b0 and b1
b0[1] = 0
b1[1] = 0


###################################################################
#### This is the part we modify (i.e. logit to modified logit) ####
loglikelihood <- function(b0, b1){
  ## Calculate p through link function
  p = exp((b0+b1*dosage))/((1 + exp(b0 + b1*dosage)))
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
  
  ## Setting prior densities
  ## I didn't have prior knowledge of beetle deaths, so I used a flat prior
  ## i.e. a normal distribution with a very large standard deviation
  current.prior.dens = dnorm(b0.prime,0,100000,log = TRUE) + dnorm(b1.prime,0,100000,log=TRUE)
  prime.prior.dens = dnorm(b0[i-1],0,100000, log = TRUE) + dnorm(b1[i-1],0,100000, log = TRUE)
  
  ## Calculating log likelihood for both the proposed and current values of b0 and b1
  current.loglik = loglikelihood(b0[i-1], b1[i-1])
  prime.loglik = loglikelihood(b0.prime, b1.prime)
  
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
  } else {
    ## Updating the chain if chain proposed value is not accepted (chain doesn't move)
    b0[i] = b0[i-1]
    b1[i] = b1[i-1]
  }
}

## Create a data frame which holds the chains for both parameters
results.vec <- data.frame(b0, b1)

layout(matrix(1:2, ncol = 1), widths = 1,
       heights = c(5,5), respect = FALSE)
par(mar=c(2.5, 4, 1, 1),oma=c(0,1,2,0))

plot(results.vec[, 1], type = "l", ylab = expression(beta[1]))
mtext("Metropolis-Hastings Traceplot")
plot(results.vec[, 2], type = "l", ylab = expression(beta[2]))



## These tend to change from sample to sample because of randomization
b0.burnin = 100000
b1.burnin = 100000

plot(results.vec[-(1:b0.burnin), 1], type = "l", ylab = expression(beta[1]))
mtext("Metropolis-Hastings Traceplot")
plot(results.vec[-(1:b1.burnin), 2], type = "l", ylab = expression(beta[2]))


## Discarding burn-in samples
b0.new <- b0[-(1:b0.burnin)]
b1.new <- b1[-(1:b1.burnin)]

# ## Checking autocorrelation using ACF
acf(b0.new, lag.max = 500)
acf(b1.new, lag.max = 500)

## From ACF on both parameters, to thin both parameters, I will take roughly every 110th sample for the b0
## and every 120th sample for b1. This is because the autocorrelation dies out after those many lags.

# ## These tend to change from sample to sample because of randomization
# b0.thin = 110 
# b1.thin = 120
# alpha.thin = 100
# 
# ## Thinning the chains
# b0.new <- b0.new[seq(1, length(b0.new), b0.thin)]
# b1.new <- b1.new[seq(1, length(b1.new), b1.thin)]
# alpha.new <- alpha.new[seq(1, length(alpha.new), alpha.thin)]

## Displaying summary stats
summary(b0.new)
summary(b1.new)

# par(mfr)
# ## Displaying thinned histograms for each parameter
# hist(b0.new, freq = FALSE, breaks = 100)
# hist(b1.new, freq = FALSE, breaks = 100)
# hist(alpha.new, freq = FALSE, breaks = 100)
