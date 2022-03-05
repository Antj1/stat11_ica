library("quantmod") # load the quantmod library
library("tseries") # load the tseries library

price = get.hist.quote(instrument = "^gspc", start = "2000-01-03",end = "2020-01-04", quote="AdjClose") # download the prices,from January 1, 2000 until today
y=diff(log(price)) # convert the prices into returns
plot(y) # plot the returns
z = coredata(y)
