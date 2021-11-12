
dat <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")


earthquakes <- ts(dat[,2], start = 1900, frequency = 1)
library(tsHMM)

data("earthquakes", package = "tsHMM")

m = 2

lambda <- c(15,25)

gamma <- matrix(c(0.9,0.1,0.1,0.9), m, m,byrow=TRUE)

model <- tsHMM$new(x = earthquakes, h = 12, m = 2, lambda = lambda, gamma = gamma)
model$fit()

forecasts <- model$forecast(xf = 0:50)

forecasts
level <- 95

# from forecast.nnetar line 84

lower <- apply(forecasts, 2, quantile, 0.5 - level/200, type = 8,
               na.rm = TRUE)
upper <- apply(forecasts, 2, quantile, 0.5 + level/200, type = 8,
               na.rm = TRUE)

upper <- apply(forecasts, 2, mean, 0.5 + level/200, type = 8,
               na.rm = TRUE)

plot(0:50, upper,type="h")



h<- 12
xf<-0:50
n   <-length(earthquakes)
par(mfrow=c(2,2),las=1)
for (i in 1:12){
  fc<-forecasts[i,]
  plot(xf,fc,type="h",main=paste("Forecast distribution for", d[n]+i),
       xlim=c(0,max(xf)),ylim=c(0,0.12),xlab="count",ylab="probability",lwd=3)
}

