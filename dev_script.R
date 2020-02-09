library(devtools)

devtools::load_all()

?HMM::HMM.pn2pw

?HMM::HMM.pw2pn

?HMM::HMM.conditional()

### A.2.1 Fitting Poisson{HMMs to the earthquakes series
dta <- read.table("http://www.hmms-for-time-series.de/second/dtaa/earthquakes.txt")
#(or set your own path)
x <- dta[,2]
d <- dta[,1]
n <- length(x)
#====================================== fit 2-state HMM
m<-2

lambda0<-c(15,25)

gamma0<-matrix(
  c(
    0.9,0.1,
    0.1,0.9
  ),m,m,byrow=TRUE)

mod2s<- HMM.mle(x,m,lambda0,gamma0,delta0 = NULL,stationary=TRUE)

delta0<-c(1,1)/2

mod2h<- HMM.mle(x,m,lambda0,gamma0,delta = delta0,stationary=FALSE)

mod2s; mod2h

#====================================== fit 3-state HMM

m<-3

lambda0<-c(10,20,30)

gamma0<-matrix(
  c(
    0.8,0.1,0.1,
    0.1,0.8,0.1,
    0.1,0.1,0.8
  ),m,m,byrow=TRUE)

mod3s<- HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)

delta0 <- c(1,1,1)/3

mod3h<- HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE)

mod3s
mod3h
#====================================== fit 4-state HMM
m <- 4

lambda0<-c(10,15,20,30)

gamma0<-matrix(
  c(
    0.85,0.05,0.05,0.05,
    0.05,0.85,0.05,0.05,
    0.05,0.05,0.85,0.05,
    0.05,0.05,0.05,0.85
  ),m,m,byrow=TRUE)

mod4s<-HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)

delta0<-c(1,1,1,1)/4

mod4h<-HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE)

mod4s; mod4h



### A.2.2 Forecast probabilities
#=== Use it for 1-step-ahead and plot the forecast distribution.
h<-1

xf<- 0:50

forecasts<-HMM.forecast(xf,h,x,mod3s)

fc<-forecasts[1,]
par(mfrow=c(1,1),las=1)
plot(xf,fc,type="h",
     main=paste("Earthquake series: forecast distribution for", d[n]+1),
     xlim=c(0,max(xf)),ylim=c(0,0.12),xlab="count",ylab="probability",lwd=3)

#=== Forecast 1-4 steps ahead and plot these.
h<-4
xf<- 0:45
forecasts<-HMM.forecast(xf,h,x,mod3s)

par(mfrow=c(2,2),las=1)
for (i in 1:4)
{
  fc<-forecasts[i,]
  plot(xf,fc,type="h",main=paste("Forecast distribution for", d[n]+i),
       xlim=c(0,max(xf)),ylim=c(0,0.12),xlab="count",ylab="probability",lwd=3)
}

#=== Compute the marginal distribution (called "dstat" below)
#    for mod3h.
#=== This is also the long-term forecast.
m<-3.

lambda<-mod3h$lambda
delta<-solve(t(diag(m)-mod3h$gamma+1),rep(1,m))
dstat<-numeric(length(xf))
for (j in 1:m) dstat <- dstat + delta[j]*dpois(xf,lambda[j])

#=== Compare the 50-year-ahead forecast with the long-term forecast.
h<-50
xf<-0:45
forecasts<-HMM.forecast(xf,h,x,mod3h)
fc<-forecasts[h,]
par(mfrow=c(1,1),las=1)
plot(xf,fc,type="h",
     main=paste("Forecast distribution for", d[n]+h),
     xlim=c(0,max(xf)),ylim=c(0,0.12),xlab="count",ylab="probability",lwd=3)
lines(xf,dstat,col="gray",lwd=3)

