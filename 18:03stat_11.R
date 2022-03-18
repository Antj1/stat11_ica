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

start_date<-as.Date("2015-01-01")
end_date<-as.Date("2021-01-01")

#Log-Returns (continuously compounded) for S&P500 & FTSE100.
priceSP = get.hist.quote(instrument="^gspc", start = start_date, end = end_date, quote="AdjClose",compression = "w")
priceFTSE = get.hist.quote(instrument="^ftse", start = start_date, end = end_date, quote="AdjClose",compression = "w")

y_individual_price = na.omit(cbind(priceSP, priceFTSE))#Prices of each asset 
y_individual = na.omit(diff(log(cbind(priceSP, priceFTSE))))#Log returns of each asset 
ggplot(data=y_individual,aes(x=Index,y=y_individual$Adjusted.priceSP))+geom_line(colour="blue")+labs(x ="t", y = "SP Log-Returns")

#Linear Correlation Matrix (Pearson's correlation coefficient).
cor(y_individual, use="pairwise.complete.obs")

#Constructing an equally weighted portfolio with prices
SP_price_init =drop(coredata(y_individual_price[1]$Adjusted.priceSP))
FTSE_price_init =drop(coredata(y_individual_price[1]$Adjusted.priceFTSE))
portfolio_weights<- c(FTSE_price_init/SP_price_init,1)

equally_weighted_price <- vector(mode = "numeric", length = dim(y_individual_price)[1])
for(i in 1:dim(y_individual_price)[1]){
  a = rowSums(y_individual_price[i]*portfolio_weights)
  equally_weighted_price[i] <- a}

dates = fortify.zoo(y_individual_price[,0])
equally_weighted_price_df = data.frame(dates,equally_weighted_price)
equally_weighted_price.xts <- xts(x=equally_weighted_price_df[,2],order.by= (equally_weighted_price_df$Index))

equally_weighted_prices = diff(y_individual_price)
net_y1 = portfolio_weights[1]*drop(coredata(equally_weighted_prices[,1]))
net_y2 = portfolio_weights[2]*drop(coredata(equally_weighted_prices[,2]))

#Individual Log-Returns (continuously compounded).
log_y1 = portfolio_weights[1]*drop(coredata(y_individual[,1]))
log_y2 = portfolio_weights[2]*drop(coredata(y_individual[,2]))

#Portfolio Log-Returns (continuously compounded).
equally_weighted_LR <- vector(mode = "numeric", length = dim(y_individual)[1])
for(i in 1:dim(y_individual)[1]){
  a = rowSums(y_individual[i])/2
  equally_weighted_LR[i] <- a}

dates = fortify.zoo(y_individual[,0])
equally_weighted_LR_df = data.frame(dates,equally_weighted_LR)
equally_weighted_LR.xts <- xts(x=equally_weighted_LR_df[,2],order.by= (equally_weighted_LR_df$Index))

#Plot of Historical Log-Returns for given Portfolio.
ggplot(data=equally_weighted_LR.xts,aes(x=Index,y=equally_weighted_LR.xts))+geom_line(colour="blue")+labs(x ="t", y = "Port Log-Returns")

#____________________________________________________________________________________

#Modelling mean (AR Model)

ret1=drop(coredata(y_individual[,1]))
ret2=drop(coredata(y_individual[,2]))

jarqueberaTest(ret1)
jarqueberaTest(ret2)

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

# Step 2: Estimation
model1=garchFit(formula=~arma(1,0)+garch(1,1),data=ret1,trace=F,cond.dist="norm")
model2=garchFit(formula=~arma(3,0)+garch(1,1),data=ret2,trace=F,cond.dist="norm")
# Note: this function uses all first three lags, not lag 3 only

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
u1<-pnorm(res1, mean=0, sd=1)[4:length(ret1)]
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
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 3)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 3)
model2@fit$ics
u2<-pnorm(res2, mean=0, sd=1)[4:length(ret2)]
hist(u2)

# Further distributional checks
#Kolmogorov-Smirnov test
KStest2<-LcKS(u2, cdf = "punif")
KStest2$p.value
#Anderson-Darling test
ADtest2<-ad.test(u2, null="punif")
ADtest2$p.value

# Misspecification Example 1: Using AIC and BIC
# Here I'm deliberately choosing the "wrong" AR(1) model instead of AR(3) model
# to show that the AIC and BIC are indeed helpful model selection criteria!
model2b=garchFit(formula=~arma(1,0)+garch(1,1),data=ret2,trace=F,cond.dist="norm")
res2b <- residuals(model2b, standardize=TRUE)
par(mfrow=c(1,2))
acf(res2b, col="green", lwd=2)
acf(res2b^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2b, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2b^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

u2b<-pnorm(res2b, mean=0, sd=1)
hist(u2b)

#Kolmogorov-Smirnov test
KStest2b<-LcKS(u2b, cdf = "punif")
KStest2b$p.value
#Anderson-Darling test
ADtest2b<-ad.test(u2b, null="punif")
ADtest2b$p.value

# Misspecification Example 2: Using Student-t instead of Normal
# Here I'm deliberately choosing the "wrong" Student-t distribution instead of 
# a Normal distribution. We know that as the df parameter approaches infinity,
# the Student-t distribution approaches the Normal distribution.
# So for the sufficiently large value of (estimated) df parameter, 
# we would not expect to see a big difference between the selection criteria.
model2c=garchFit(formula=~arma(3,0)+garch(1,1),data=ret2,trace=F, cond.dist="std")
res2c <- residuals(model2c, standardize=TRUE)
par(mfrow=c(1,2))
acf(res2c, col="green", lwd=2)
acf(res2c^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2c, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2c^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
coef(model2c)

u2c<-pnorm(res2c, mean=0, sd=1)
hist(u2c)

#Kolmogorov-Smirnov test
KStest2c<-LcKS(u2c, cdf = "punif")
KStest2c$p.value
#Anderson-Darling test
ADtest2c<-ad.test(u2c, null="punif")
ADtest2c$p.value

ic=rbind(model2@fit$ics,model2b@fit$ics,model2c@fit$ics)
rownames(ic)<-c("AR(3)","AR(1)","Student-t")
ic

# Let's compare p-values for the KS and AD tests for uniformity
dtests=rbind(c(KStest2$p.value,ADtest2$p.value),c(KStest2b$p.value,ADtest2b$p.value),c(KStest2c$p.value,ADtest2c$p.value))
rownames(dtests)<-c("AR(3)","AR(1)","Student-t")
dtests

# Copula modelling
model=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
model

# Value-at-Risk uisng Monte Carlo simulation
N=10000
u_sim=BiCopSim(N, family=model$family, model$par,  model$par2)
# Here we are assuming marginal models are N(0,1), completely ignoring our knowledge of real DGP
res1_sim=qnorm(u_sim[,1], mean = 10, sd = 4) 
res2_sim=qnorm(u_sim[,2], mean = 8, sd = 7) 

# However, "res1_sim" and "res2_sim" are i.i.d.
# So, the next step is to re-introduce  autocorrelation and GARCH effects observed in data

# This section has been omitted as this will be part of your ICA group assignment
#####################################################################################




#####################################################################################
portsim <- matrix(0, nrow = N, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

#Using log_y1 e.t.c is calling past returns, we need to simulate future returns instead
portsim=log(1+((exp(log_y1)-1)+(exp(log_y2)-1))*(1/2))
portsim
varsim=last(equally_weighted_price,1)*quantile(portsim,c(0.01,0.05))
varsim

