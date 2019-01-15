data = read.csv("gator2.txt", header = T, sep = "")
with(data, table(choice, sex))

data$choice <- relevel(data$choice, ref = "F")
data$choice = as.numeric(data$choice)

log(data$length)

# Code 5.16 Logistic model using R
library(MCMCpack)
# ?MCMCmnl
posterior.mcmc <- MCMCmnl(choice ~ length,
               burnin = 5000,
               mcmc = 10000,
               data = data)




par(mar=c(3,5,2,1))
caterplot(myL1)
densplot(myL1[,1], main = "Posterior Density Plot of Intercept\n Estimate for Invertebrate Category")
densplot(myL1[,2], main = "Posterior Density Plot of Intercept\n Estimate for Other Category")
densplot(myL1[,3], main = "Posterior Density Plot of Beta\n Estimate for Invertebrate Category")
densplot(myL1[,4], main = "Posterior Density Plot of Beta\n Estimate for Other Category")
traceplot(myL1[,1],main = "Trace Plot of Intercept Estimate\n for Invertebrate Category")
traceplot(myL1[,2],main = "Trace Plot of Intercept Estimate\n for Other Category")
traceplot(myL1[,3],main = "Trace Plot of Beta Estimate\n for Invertebrate Category")
traceplot(myL1[,4],main = "Trace Plot of Intercept Estimate\n for Other Category")

### summary statistics
summaryout = cbind(summary(myL1)$statistics[,c(1,2)],summary(myL1)$quantiles)
# library(xtable)
xtable(summaryout)

pi_11 = function(x1){
  num = exp(4.6596082 - 2.7579199*x1)
  den = 1 + exp(4.6596082 - 2.7579199*x1) + exp(-0.9621408 - 0.1320392*x1)
  return(num/den)
}
pi_21 = function(x1){
  num = exp(-0.9621408 - 0.1320392*x1)
  den = 1 + exp(4.6596082 - 2.7579199*x1) + exp(-0.9621408 - 0.1320392*x1)
  return(num/den)
}
pi_31 = function(x1){
  num = 1
  den = 1 + exp(4.6596082 - 2.7579199*x1) + exp(-0.9621408 - 0.1320392*x1)
  return(num/den)
}

x1grid = seq(0,4,0.01)
plot(x1grid, pi_11(x1grid), col = "white", xlab = "exp(length)", ylim = c(0,1),
     ylab = expression(hat(pi[i])),main = "Predicted Probabilities")
lines(x1grid, pi_11(x1grid), lwd=1, col = "red")
lines(x1grid, pi_21(x1grid), col = "blue")
lines(x1grid, pi_31(x1grid), col = "green")
legend("left",legend=c('Invertebrates','Other','Fish'),
       text.col=c("red", "blue", "green"),cex=1)


######################## Model 2 ######################################
myL <- MCMCmnl(choice ~ sex + length,
                 burnin = 10000,
                 mcmc = 100000,
                 data = data)
summaryout2 = cbind(summary(myL)$statistics[,c(1,2)],summary(myL)$quantiles)
summaryout2
# library(xtable)
xtable(summaryout2)

par(mar=c(3,5,3,1))
caterplot(myL)
mtext("Plot of Credible Intervals for\n Parameters from MCMC Simulation", side = 3, font =2)

?caterplot

summary(myL)
plot.ts(myL[,1])
plot(density(myL[,1]))
library(coda)



pi_1 = function(x1,x2){
  num = exp(6.399 - 1.33*x1 - 3.29*x2)
  den = 1 + exp(6.399 - 1.33*x1 - 3.29*x2) + exp(-1.1581 + 0.2019*x1 - 0.1233*x2)
  return(num/den)
}
pi_2 = function(x1,x2){
  num = exp(-1.1581 + 0.2019*x1 - 0.1233*x2)
  den = 1 + exp(6.399 - 1.33*x1 - 3.29*x2) + exp(-1.1581 + 0.2019*x1 - 0.1233*x2)
  return(num/den)
}
pi_3 = function(x1,x2){
  num = 1
  den = 1 + exp(6.399 - 1.33*x1 - 3.29*x2) + exp(-1.1581 + 0.2019*x1 - 0.1233*x2)
  return(num/den)
}

##### Males ######
x1grid = seq(0,4,0.01)
plot(x1grid, logit1(1, x1grid), col = "white", xlab = "exp(length)", ylim = c(0,1),
     ylab = expression(hat(pi[i])),main = "Predicted probability Male Alligators")
lines(x1grid, pi_1(1, x1grid), lwd=1, col = "red")
lines(x1grid, pi_2(1, x1grid), col = "blue")
lines(x1grid, pi_3(1, x1grid), col = "green")
legend("right",legend=c('Invertebrates','Other','Fish'),
       text.col=c("red", "blue", "green"),cex=1)

x2grid = seq(0,4,0.01)
plot(x2grid, logit1(0,x2grid), col = "white", xlab = "exp(length)", ylim = c(0,1),
     ylab = expression(hat(pi[i])),main = "Predicted probability Female Alligators")
lines(x2grid, pi_1(0,x2grid), lwd=1, col = "orange")
lines(x2grid, pi_2(0, x2grid), col = "purple")
lines(x2grid, pi_3(0, x2grid), col = "cyan")

legend("right",legend=c('Invertebrates','Other','Fish'),
       text.col=c("orange", "purple", "cyan"),cex=1)

##### plots overlayed #####
x1grid = seq(0,4,0.01)
plot(x1grid, pi_1(1, x1grid), col = "white", xlab = "exp(length)", ylim = c(0,1),
     ylab = expression(hat(pi[i])),main = "Predicted probability Overlayed\n (to see the difference)")
lines(x1grid, pi_1(1, x1grid), lwd=1, col = "red")
lines(x1grid, pi_2(1, x1grid), col = "blue")
lines(x1grid, pi_3(1, x1grid), col = "green")
lines(x1grid, pi_1(0,x1grid), lwd=1, col = "orange")
lines(x1grid, pi_2(0, x1grid), col = "purple")
lines(x1grid, pi_3(0, x1grid), col = "cyan")


##### Point and interval estimate plots model 1 ######
# Thus, we have a distribution of model predictions for each x
predict.p.mcmc <- function(x) 1 / (1 + exp(-posterior.mcmc[,c(1,2,3,4)] %*% rbind(1,x)))
interval.p.mcmc <- function(x, low, high) apply(posterior.mcmc[,c(1,2,3,4)](x), 2,
                                                function(x) quantile(x, c(low, high)))

predict.y.mcmc <- function(x) posterior.mcmc %*% rbind(1,x)
interval.y.mcmc <- function(x, low, high) apply(predict.y.mcmc(x), 2,
                                                function(x) quantile(x, c(low, high)))
x_test = seq(0,4, length.out = 200)
interval.p.mcmc_95 <- interval.p.mcmc(x_test, 0.025, 0.975)
lines(x_test, interval.p.mcmc_95[1,], col = 'red')
lines(x_test, interval.p.mcmc_95[2,], col = 'red')



summary(myL1)














