library("tseries")
library("zoo")

library(PerformanceAnalytics)

#Log (continuously compounded) Returns for S&P500 and FTSE100
priceSP = get.hist.quote(instrument="^gspc", start = "2010-03-09", end = "2020-12-15", quote="AdjClose")
y_sp=diff(log(priceSP))
plot(y_sp)
y_sp=coredata(y_sp)

priceFTSE = get.hist.quote(instrument="^ftse", start = "2010-03-09", end = "2020-12-15", quote="AdjClose")
y_ftse=diff(log(priceFTSE))
plot(y_ftse)
y_ftse=coredata(y_ftse)

#Stock prices, and indices correlation matrix
priceGE = get.hist.quote(instrument="ge", start = "2010-03-09", end = "2020-12-15", quote="AdjClose")
priceMSFT = get.hist.quote(instrument="msft", start = "2010-03-09", end = "2020-12-15", quote="AdjClose")

p <- cbind(priceSP, priceFTSE, priceGE, priceMSFT)
y_individual = diff(log(p))#Log returns of each asset 
y_individual = na.omit(y_individual)
cor(y_individual, use="pairwise.complete.obs")

#Now constructing an equally weighted portfolio with log returns


equal_weights <- vector(mode = "numeric", length = dim(y_individual)[1])

dates = fortify.zoo(y_individual[,0])


for(i in 1:dim(y_individual)[1]){
  a = rowSums(y_individual[i])/4
  equal_weights[i] <- a}

plot(equalp)
plot(dates,equal_weights,type = "l")

