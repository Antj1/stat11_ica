library("tseries")
library("zoo")

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
p=cbind(priceSP, priceFTSE, priceGE, priceMSFT)
y_port = diff(log(p))
cor(y_port, use="pairwise.complete.obs")
