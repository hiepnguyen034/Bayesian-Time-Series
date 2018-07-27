require(msm)
library(MCMCpack)
library(forecast)
library(stargazer)
library(tseries)
set.seed(2020)

x <- arima.sim(model = list(order = c(1, 0, 0), ar = .3), n = 100) 
arma(x,p=1,q=0)
pacf(x)

yt       <- x
ytmin1   <- yt[c(1:(length(yt)-1))]
yt       <- yt[c(2:length(yt))]
prod     <- yt%*%ytmin1
prod     <- prod[1]
ytsq     <- sum(ytmin1^2)
n        <- length(yt)

B		       <- 1000
phi		     <- vector("numeric", B)
mu	       <- vector("numeric", B)
sig2eps	   <- vector("numeric", B)
sig2phi    <- vector("numeric", B)

phi[1]            <- mean(yt)
mu[1]             <- 0
sig2eps[1]        <- var(yt)
sig2phi[1]        <- 1

for(t in 2:B){
  ### sample phi ###
  phi[t] <- rnorm(1, (prod/sig2eps[t-1]+mu[t-1]/sig2phi[t-1])*(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1,(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1)
  ### sample mu ###
  mu[t] <- rtnorm(1,phi[t-1],sig2phi[t-1],-1,1)
  ### sample sig2eps ###
  sig2eps[t] <- rinvgamma(1,n/2,(1/2)*sum((yt-phi[t-1]*ytmin1)^2))
  ### sample sig2phi ###
  sig2phi[t] <- rinvgamma(1,1/2,(1/2)*(phi[t-1]-mu[t-1])^2)
}

par(mfrow=c(2,2))

plot.ts(phi)
plot.ts(mu)
plot.ts(sig2eps)
plot.ts(sig2phi)

##############

library(astsa)
require(astsa)

data(gtemp)

plot.ts(gtemp)

diff.gtemp <- diff(gtemp)

############## 


yt       <- diff.gtemp 
ytmin1   <- yt[c(1:(length(yt)-1))]
yt       <- yt[c(2:length(yt))]
prod     <- yt%*%ytmin1
prod     <- prod[1]
ytsq     <- sum(ytmin1^2)
n        <- length(yt)

B		       <- 20000
phi		     <- vector("numeric", B)
mu	       <- vector("numeric", B)
sig2eps	   <- vector("numeric", B)
sig2phi    <- vector("numeric", B)

phi[1]            <- mean(yt)
mu[1]             <- 0
sig2eps[1]        <- var(yt)
sig2phi[1]        <- 1

for(t in 2:B){
  ### sample phi ###
  phi[t] <- rnorm(1, (prod/sig2eps[t-1]+mu[t-1]/sig2phi[t-1])*(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1,(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1)
  ### sample mu ###
  mu[t] <- rtnorm(1,phi[t-1],sig2phi[t-1],-1,1)
  ### sample sig2eps ###
  sig2eps[t] <- rinvgamma(1,n/2,(1/2)*sum((yt-phi[t-1]*ytmin1)^2))
  ### sample sig2phi ###
  sig2phi[t] <- rinvgamma(1,1/2,(1/2)*(phi[t-1]-mu[t-1])^2)
}
warnings()

par(mfrow=c(2,2))

plot.ts(phi)
plot.ts(mu)
plot.ts(sig2eps)
plot.ts(sig2phi)

pacf(phi[10000:20000])
pacf(mu[10000:20000])
pacf(sig2eps[10000:20000])
pacf(sig2phi[10000:20000])

plot(density(phi[10000:20000]))
plot(density(mu[10000:20000]))
plot(density(sig2eps[10000:20000]))
plot(density(sig2phi[10000:20000]))


geweke.diag(phi)
geweke.diag(mu)
geweke.diag(sig2eps)
geweke.diag(sig2phi)

gelman.diag(phi)
gelman.diag(mu)
gelman.diag(sig2eps)
gelman.diag(sig2phi)


most.thin <- rep(c(rep(F, 19), T), 1000)

summary(most.thin)
length(phi)

pacf(phi[most.thin])
pacf(mu[most.thin])
pacf(sig2eps[most.thin])
pacf(sig2phi[most.thin])

tab <- t(cbind(quantile(phi[most.thin],     probs=c(.025, .5, .975)), 
               quantile(mu[most.thin],      probs=c(.025, .5, .975)), 
               quantile(sig2eps[most.thin], probs=c(.025, .5, .975)), 
               quantile(sig2phi[most.thin], probs=c(.025, .5, .975))))

stargazer(tab) 

model.frequentist <- arma(diff.gtemp, order = c(1,0))

###################

# Deseasonalized Brazil GDP

###################


brazil2 <- read.csv("./brazil_predictors.txt", header = T)

head(brazil2)

brazil.GDP <- ts(brazil2$GDP, frequency=12)

decomposed.brazil <- decompose(brazil.GDP)

decomposed.brazil

brazil.deseasonalized <- na.omit(decomposed.brazil$trend + decomposed.brazil$random)

plot.ts(brazil.deseasonalized) 

brazil.deseasonalized.diff <- vector("numeric", length(brazil.deseasonalized))

for(i in 2:length(brazil.deseasonalized)){
  
  brazil.deseasonalized.diff[i] <- brazil.deseasonalized[i] - lag(brazil.deseasonalized[i-1], 1)
  
}



yt       <- brazil.deseasonalized.diff 
ytmin1   <- yt[c(1:(length(yt)-1))]
yt       <- yt[c(2:length(yt))]
prod     <- yt%*%ytmin1
prod     <- prod[1]
ytsq     <- sum(ytmin1^2)
n        <- length(yt)

B		       <- 20000
phi		     <- vector("numeric", B)
mu	       <- vector("numeric", B)
sig2eps	   <- vector("numeric", B)
sig2phi    <- vector("numeric", B)

phi[1]            <- mean(yt)
mu[1]             <- 0
sig2eps[1]        <- var(yt)
sig2phi[1]        <- 1

for(t in 2:B){
  ### sample phi ###
  phi[t] <- rnorm(1, (prod/sig2eps[t-1]+mu[t-1]/sig2phi[t-1])*(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1,(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1)
  ### sample mu ###
  mu[t] <- rtnorm(1,phi[t-1],sig2phi[t-1],-1,1)
  ### sample sig2eps ###
  sig2eps[t] <- rinvgamma(1,n/2,(1/2)*sum((yt-phi[t-1]*ytmin1)^2))
  ### sample sig2phi ###
  sig2phi[t] <- rinvgamma(1,1/2,(1/2)*(phi[t-1]-mu[t-1])^2)
}

par(mfrow=c(2,2))

plot.ts(phi)
plot.ts(mu)
plot.ts(sig2eps)
plot.ts(sig2phi)

pacf(phi[10000:20000])
pacf(mu[10000:20000])
pacf(sig2eps[10000:20000])
pacf(sig2phi[10000:20000])

plot(density(phi[10000:20000]))
plot(density(mu[10000:20000]))
plot(density(sig2eps[10000:20000]))
plot(density(sig2phi[10000:20000]))

geweke.diag(phi)
geweke.diag(mu)
geweke.diag(sig2eps)
geweke.diag(sig2phi)

# Gelman diagnostics for phi

first.chain  <- rep(c(T,F,F,F), B/4)
second.chain <- rep(c(F,T,F,F), B/4)
third.chain  <- rep(c(F,F,T,F), B/4)
fourth.chain <- rep(c(F,F,F,T), B/4)

chain1	<- mcmc(phi[first.chain])
chain2  <- mcmc(phi[second.chain])
chain3  <- mcmc(phi[third.chain])
chain4  <- mcmc(phi[fourth.chain])

allChains <- mcmc.list(list(chain1, chain2, chain3, chain4)) 

gelman.diag(allChains)

# Gelman diagnostics for mu


chain1	<- mcmc(mu[first.chain])
chain2  <- mcmc(mu[second.chain])
chain3  <- mcmc(mu[third.chain])
chain4  <- mcmc(mu[fourth.chain])

allChains <- mcmc.list(list(chain1, chain2, chain3, chain4)) 

gelman.diag(allChains)


most.thin <- rep(c(rep(F, 19), T), 1000)

summary(most.thin)
length(phi)

pacf(phi[most.thin])
pacf(mu[most.thin])
pacf(sig2eps[most.thin])
pacf(sig2phi[most.thin])

tab <- t(cbind(quantile(phi[most.thin],     probs=c(.025, .5, .975)), 
               quantile(mu[most.thin],      probs=c(.025, .5, .975)), 
               quantile(sig2eps[most.thin], probs=c(.025, .5, .975)), 
               quantile(sig2phi[most.thin], probs=c(.025, .5, .975))))

stargazer(tab) 

model.frequentist.2 <- arma(brazil.deseasonalized.diff, order = c(1,0))


##############

# Brazil GDP -- nonstationary

############## 

brazil <- read.table("./brazilgdp.txt")

plot.ts(brazil$V1)

yt       <- brazil$V1 
ytmin1   <- yt[c(1:(length(yt)-1))]
yt       <- yt[c(2:length(yt))]
prod     <- yt%*%ytmin1
prod     <- prod[1]
ytsq     <- sum(ytmin1^2)
n        <- length(yt)

B		       <- 10000
phi		     <- vector("numeric", B)
mu	       <- vector("numeric", B)
sig2eps	   <- vector("numeric", B)
sig2phi    <- vector("numeric", B)

phi[1]            <- 1
mu[1]             <- 0
sig2eps[1]        <- 100
sig2phi[1]        <- 1

for(t in 2:B){
  ### sample phi ###
  phi[t] <- rnorm(1, (prod/sig2eps[t-1]+mu[t-1]/sig2phi[t-1])*(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1,(ytsq/sig2eps[t-1]+1/sig2phi[t-1])^-1)
  ### sample mu ###
  mu[t] <- rtnorm(1,phi[t-1],sig2phi[t-1], -100, 100)
  ### sample sig2eps ###
  sig2eps[t] <- rinvgamma(1,n/2,(1/2)*sum((yt-phi[t-1]*ytmin1)^2))
  ### sample sig2phi ###
  sig2phi[t] <- rinvgamma(1,1/2,(1/2)*(phi[t-1]-mu[t-1])^2)
}

par(mfrow=c(2,2))

plot.ts(phi)
plot.ts(mu)
plot.ts(sig2eps)
plot.ts(sig2phi)

pacf(phi[B/2:B])
pacf(mu[B/2:B])
pacf(sig2eps[B/2:B])
pacf(sig2phi[B/2:B])

plot(density(phi[B/2:B]))
plot(density(mu[B/2:B]))
plot(density(sig2eps[B/2:B]))
plot(density(sig2phi[B/2:B]))

geweke.diag(phi)
geweke.diag(mu)
geweke.diag(sig2eps)
geweke.diag(sig2phi)

# Gelman diagnostics for phi

first.chain  <- rep(c(T,F,F,F), B/4)
second.chain <- rep(c(F,T,F,F), B/4)
third.chain  <- rep(c(F,F,T,F), B/4)
fourth.chain <- rep(c(F,F,F,T), B/4)

chain1	<- mcmc(phi[first.chain])
chain2  <- mcmc(phi[second.chain])
chain3  <- mcmc(phi[third.chain])
chain4  <- mcmc(phi[fourth.chain])

allChains <- mcmc.list(list(chain1, chain2, chain3, chain4)) 

gelman.diag(allChains)


# Gelman diagnostics for mu

chain1	<- mcmc(mu[first.chain])
chain2  <- mcmc(mu[second.chain])
chain3  <- mcmc(mu[third.chain])
chain4  <- mcmc(mu[fourth.chain])

allChains <- mcmc.list(list(chain1, chain2, chain3, chain4)) 

gelman.diag(allChains)

most.thin <- rep(c(rep(F, 19), T), B/20)

summary(most.thin)
length(phi)

pacf(phi[most.thin])
pacf(mu[most.thin])
pacf(sig2eps[most.thin])
pacf(sig2phi[most.thin])

tab <- t(cbind(quantile(phi[most.thin],     probs=c(.025, .5, .975)), 
               quantile(mu[most.thin],      probs=c(.025, .5, .975)), 
               quantile(sig2eps[most.thin], probs=c(.025, .5, .975)), 
               quantile(sig2phi[most.thin], probs=c(.025, .5, .975))))

stargazer(tab) 


model.frequentist.2 <- arma(brazil$V1, order = c(1,0))

summary(model.frequentist.2)
