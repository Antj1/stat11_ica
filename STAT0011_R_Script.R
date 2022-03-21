###########
# TASK 1 #
###########

data = read.csv("/Users/harrisonknowles/Desktop/FINISHED_STAT0011/group27.csv")
x = data[,1]
i = 1:299
prob_tau = 1/298

integranda <- function(sigsq,tau) {
  u<-1
  for(i in 1:tau){
    u = u * 1/(sqrt(2*pi*sigsq)) * exp(-1/2*(x[i]/sqrt(sigsq))^2)
  }
  return(u)
}

integranda2 <- function(sigsq,tau) {
  u<-1
  for(i in (tau:298)){
    u = u * 1/(sqrt(2*pi*sigsq)) * exp(-1/2*(x[i]/sqrt(sigsq))^2)
  }
  return(u)
}

integrandb <- function(sigsq){
  v = (1/(sigsq))^2 * exp(-1/sigsq)
  return(v)
}

integrand1 <- function(sigsq,tau){
  int = integranda(sigsq,tau)*integrandb(sigsq)
  return(int)
}

integrand2 <- function(sigsq,tau){
  int = integranda2(sigsq,tau)*integrandb(sigsq)
  return(int)
}

probs<- rep(0,298)
j<-1
while(j < length(probs)){
  probs[j]=(integrate(integrand1,0,Inf,tau=j)$value)*(integrate(integrand2,0,Inf,tau=j)$value)
  j = j+1
}
probs = (probs*prob_tau) / sum(probs)
which.max(probs)

plot(i,x)
abline(v=145, col="blue")
text(x = 120, y = -6, label = "Tau=145", srt = 0,cex=1)

###########
# TASK 2 #
###########

library("tseries")
library("zoo")
library(gridExtra)
library(ggTimeSeries)
library(xts)
library(pracma)
library(VineCopula)
library(fGarch)
library(KScorrect)
library(stats)
library(ADGofTest)
rm(list=ls(all=TRUE))
dev.off()

start_date<-as.Date("2009-10-01")
end_date<-as.Date("2016-02-02")

#Log-Returns (continuously compounded) for S&P500 & FTSE100.
priceSP = get.hist.quote(instrument="^gspc", start = start_date, end = end_date, quote="AdjClose",compression = "w")
priceFTSE = get.hist.quote(instrument="^ftse", start = start_date, end = end_date, quote="AdjClose",compression = "w")
y_individual_price = na.omit(cbind(priceSP, priceFTSE))#Prices of each asset 
y_individual = na.omit(diff(log(cbind(priceSP, priceFTSE))))#Log returns of each asset 

ggplot(data=y_individual,aes(x=Index,y=Adjusted.priceSP))+geom_line(colour="blue")+labs(x ="t", y = "SP Log-Returns")
ggplot(data=y_individual,aes(x=Index,y=Adjusted.priceFTSE))+geom_line(colour="blue")+labs(x ="t", y = "FTSE Log-Returns")

#Constructing asset weights for the equally-weighted portfolio.
SP_price_init =drop(coredata(y_individual_price[1]$Adjusted.priceSP))
FTSE_price_init =drop(coredata(y_individual_price[1]$Adjusted.priceFTSE))
portfolio_weights<- c(FTSE_price_init/SP_price_init,1)
weightsum = sum(portfolio_weights)
portfolio_weights[1] = portfolio_weights[1]/weightsum 
portfolio_weights[2] = portfolio_weights[2]/weightsum 

#Individual Historical Log-Returns.
ret1 = portfolio_weights[1]*drop(coredata(y_individual[,1])) #Log-Returns of S&P500
ret2 = portfolio_weights[2]*drop(coredata(y_individual[,2])) #Log-Returns of FTSE100

#Historical Portfolio Log-Returns.
equally_weighted_LR <- vector(mode = "numeric", length = dim(y_individual)[1])
for(i in 1:dim(y_individual)[1]){
  a = rowSums(y_individual[i])/2
  equally_weighted_LR[i] <- a}
dates = fortify.zoo(y_individual[,0])
equally_weighted_LR_df = data.frame(dates,equally_weighted_LR)
equally_weighted_LR.xts <- xts(x=equally_weighted_LR_df[,2],order.by= (equally_weighted_LR_df$Index))

#Plot of Historical Log-Returns for given Portfolio.
ggplot(data=equally_weighted_LR.xts,aes(x=Index,y=equally_weighted_LR))+geom_line(colour="blue")+labs(x ="t", y = "Port Log-Returns")

#____________________________________________________________________________________

jarqueberaTest(ret1) #test of whether sample data have the skewness and kurtosis matching a normal distribution.
jarqueberaTest(ret2)
#Small p-value implies reject Null-Hypothesis i.e. these data do not come from a Normal Distribution.

# Building AR models: the Box - Jenkins approach
# Step 1: Identification

# returns 1
par(mfrow=c(2,2))
acf(ret1, col="green", lwd=2)
pacf(ret1, col="green", lwd=2)
acf(ret1^2, col="red", lwd=2)
par(mfrow=c(1,1))

# returns 2
par(mfrow=c(2,2))
acf(ret2, col="green", lwd=2)
pacf(ret2, col="green", lwd=2)
acf(ret2^2, col="red", lwd=2)
par(mfrow=c(1,1))

# Step 2: Estimation (Variations & Results have been recorded in report table)
model1=garchFit(formula=~arma(1,0)+garch(1,1),data=ret1,trace=F,cond.dist="ged")
model2=garchFit(formula=~arma(1,0)+garch(1,1),data=ret2,trace=F,cond.dist="std")
#Lags determined from acf plots, garch parameters to be selected by AIC minimization.
#Experiment with cond.dist options, what provides best fit???

# Step 3: Model checking
# returns 1
res1 <- residuals(model1, standardize=TRUE)
par(mfrow=c(2,1))
acf(res1, col="green", lwd=2)
acf(res1^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
model1@fit$ics
coef(model1)
#u1<-pnorm(res1, mean=0, sd=1)[2:length(ret1)]
#u1<-pstd(res1, mean=0, sd=1,nu=coef(model1)[6])[2:length(ret1)]
#u1<-psstd(res1, mean=0, sd=1,nu=coef(model1)[7])[2:length(ret1)]
u1<-pged(res1, mean=0, sd=1,nu=coef(model1)[6])[2:length(ret1)]
hist(u1)

# Further distributional checks
#Kolmogorov-Smirnov test
KStest1<-LcKS(u1, cdf = "punif")
KStest1$p.value
#Anderson-Darling test
ADtest1<-ad.test(u1, null="punif")
ADtest1$p.value

# returns 2
res2 <- residuals(model2, standardize=TRUE)
par(mfrow=c(2,1))
acf(res2, col="green", lwd=2)
acf(res2^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
model2@fit$ics
coef(model2)
#u2<-pnorm(res2, mean=0, sd=1)[2:length(ret2)]
u2<-pstd(res2, mean=0, sd=1,nu=coef(model2)[6])[2:length(ret2)]
#u2<-psstd(res2, mean=0, sd=1,nu=coef(model2)[7])[2:length(ret2)]
#u2<-pged(res2, mean=0, sd=1,nu=coef(model2)[6])[2:length(ret2)]
hist(u2)

# Further distributional checks
#Kolmogorov-Smirnov test
KStest2<-LcKS(u2, cdf = "punif")
KStest2$p.value
#Anderson-Darling test
ADtest2<-ad.test(u2, null="punif")
ADtest2$p.value

#____________________
#TESTING RESULTS IN REPORT TABLE
#____________________

model1=garchFit(formula=~arma(1,0)+garch(1,1),data=ret1,trace=F,cond.dist="ged")
model2=garchFit(formula=~arma(1,0)+garch(1,1),data=ret2,trace=F,cond.dist="std")
u1<-pged(res1, mean=0, sd=1,nu=coef(model1)[6])[2:length(ret1)]
u2<-pstd(res2, mean=0, sd=1,nu=coef(model2)[6])[2:length(ret2)]
# Copula modelling
model=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
model

# Value-at-Risk uisng Monte Carlo simulation
N=10000
u_sim=BiCopSim(N, family=model$family, model$par,  model$par2)
# Here we are assuming marginal models are N(0,1), completely ignoring our knowledge of real DGP
res1_sim=qged(u_sim[,1], mean = 1, sd = 1) 
res2_sim=qstd(u_sim[,2], mean = 1, sd = 1) 

uu1 = u_sim[,1]
uu2 = u_sim[,2]

hist(u_sim[,1])
hist(u_sim[,2])
hist(res1_sim)
hist(res2_sim)
plot(u_sim)

# However, "res1_sim" and "res2_sim" are i.i.d.
acf(res1_sim) #We see no significant lags because this is a simulated i.i.d process
acf(res2_sim) #We see no significant lags because this is a simulated i.i.d process

# Re-introduction of autocorrelation and GARCH effects observed in historical data
#Simulate sigma_t N times from GARCH
OMEGA_1 = coef(model1)[3]
MU_1 = coef(model1)[1]
AR_1 = coef(model1)[2]
ALPH_1 = coef(model1)[4]
BETA_1 = coef(model1)[5]

Meansim1 <- matrix(0, nrow = N, ncol = 1)
Meansim1 [1] = MU_1 + AR_1*ret1[331]
for(i in seq(2, N, 1)) {
  Meansim1[i] = MU_1 + AR_1*Meansim1[i-1]
}

Variancesim1 <- matrix(0, nrow = N, ncol = 1)
Variancesim1 [1] = OMEGA_1 + ALPH_1*(model1@h.t[331])*(ret1[331])^2 + BETA_1*model1@h.t[331]
for(i in seq(2, N, 1)) {
  Variancesim1[i] = OMEGA_1 + ALPH_1*(uu1[i-1])^2 + BETA_1*Variancesim1[i-1]
}


Returnsim1 <- matrix(0, nrow = N, ncol = 1)
for(i in seq(1, N, 1)) {
  Returnsim1[i] = Meansim1[i]+Variancesim1[i]*res1_sim[i]
}

OMEGA_2 = coef(model2)[3]
MU_2 = coef(model2)[1]
AR_2 = coef(model2)[2]
ALPH_2 = coef(model2)[4]
BETA_2 = coef(model2)[5]

Meansim2 <- matrix(0, nrow = N, ncol = 1)
Meansim2 [1] = MU_2 + AR_2*ret2[331]
for(i in seq(2, N, 1)) {
  Meansim2[i] = MU_2 + AR_2*Meansim2[i-1]
}

Variancesim2 <- matrix(0, nrow = N, ncol = 1)
Variancesim2 [1] = OMEGA_2 + ALPH_2*(model2@h.t[331])*(ret2[331])^2 + BETA_2*model2@h.t[331]
for(i in seq(2, N, 1)) {
  Variancesim2[i] = OMEGA_2 + ALPH_2*(uu2[i-1])^2 + BETA_2*Variancesim2[i-1]
}

Returnsim2 <- matrix(0, nrow = N, ncol = 1)
for(i in seq(1, N, 1)) {
  Returnsim2[i] = Meansim2[i]+Variancesim2[i]*res2_sim[i]
}


portsim <- matrix(0, nrow = N, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)
portsim=log(1+(exp(Returnsim1)-1)*(portfolio_weights[1])+(exp(Returnsim2)-1)*(portfolio_weights[2]))
varsim= abs(quantile(portsim,c(0.01,0.05)))
varsim

