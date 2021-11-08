#' @title Forecasting Using Hidden Markov Models
#' @description
#' Forecasting Using Hidden Markov Models
#' @examples
#' \dontrun{
#'
#' library(tsHMM)
#'
#' data("earthquakes", package = "tsHMM")
#'
#' m = 2
#'
#' lambda <- c(15,25)
#'
#' gamma <- matrix(c(0.9,0.1,0.1,0.9), m, m,byrow=TRUE)
#'
#' model <- tsHMM$new(x = earthquakes, h = 12, m = 2, lambda = lambda, gamma = gamma)
#'
#' model$fit()
#'
#' fc <- model$forecast(xf = 0:50)
#'
#' fc
#'
#' }
#' @author Resul Akay
#'
#' @references
#' Zucchini, W., MacDonald, I.L., & Langrock, R. (2016).
#' Hidden Markov Models for Time Series: An Introduction Using R (2nd ed.).
#' Chapman and Hall/CRC. https://doi.org/10.1201/b20790
#' @source \url{http://www.hmms-for-time-series.de/second/appendix-a/ZMLcode.txt}
#'
#' @importFrom R6 R6Class
#' @export
tsHMM <- R6::R6Class(
  classname = "tsHMM",
  public = list(
    #' @description Initialize a HMM
    #' @param x A univariate time series
    #' @param h Forecast horizon
    #' @param m Number of stats
    #' @param lambda lambda
    #' @param gamma gamma
    initialize = function(x, h = frequency(x), m, lambda, gamma){

      if(class(x)!= "ts"){
        stop("x must be a time series")
      }

      private$x <- x
      private$d <- time(x)
      private$n <- length(x)
      private$h <- h
      private$m <- m
      private$lambda <- lambda
      private$gamma <- gamma
    },
    #' @description Fit a HMM to a univariate timeseries
    #' @param delta delta
    #' @param stationary Bool. If stationary = TRUE data.
    fit = function(delta = NULL, stationary = TRUE){
      private$stationary <- stationary
      private$delta <- delta
      private$fitted_model <- private$mle(delta = delta)
      return(invisible(NULL))
    },
    #' @description Forecast fitted HMM
    #' @param xf xf
    forecast = function(xf = NULL){

      forecasts <- private$forecast_(xf)
      return(forecasts)
    }
  ),
  private = list(
    x = NULL,
    d = NULL,
    n = NULL,
    h = NULL,
    m = NULL,
    lambda = NULL,
    gamma = NULL,
    stationary = NULL,
    fitted_model = NULL,
    delta = NULL,
    # Transforming natural parameters to working
    pn2pw = function(m, lambda, gamma, delta = NULL, stationary = TRUE){
      tlambda <- log(lambda)
      if(m==1) return(tlambda)
      foo <- log(gamma/diag(gamma))
      tgamma <- as.vector(foo[!diag(m)])
      if(stationary){
        tdelta  <- NULL
      }
      else {
        tdelta <- log(delta[-1]/delta[1])
      }
      parvect <- c(tlambda,tgamma,tdelta)
      return(parvect)
    },
    # Transforming working parameters to natural\
    pw2pn = function(m,parvect,stationary=TRUE){
      lambda <- exp(parvect[1:m])
      gamma <- diag(m)
      if (m==1) {
        return(list(lambda = lambda, gamma = gamma, delta = 1))
      }
      gamma[!gamma] <- exp(parvect[(m+1):(m*m)])
      gamma <- gamma/apply(gamma,1,sum)
      if(stationary){
        delta <- solve(t(diag(m)-gamma+1),rep(1,m))
      }
      else {
        foo<-c(1,exp(parvect[(m*m+1):(m*m+m-1)]))
        delta<-foo/sum(foo)
      }
      return(list(lambda = lambda, gamma = gamma, delta = delta))
    },
    # Computing minus the log-likelihood from the working parameters
    mllk = function(parvect,x,m,stationary=TRUE,...){
      if(m==1) return(-sum(dpois(x,exp(parvect),log=TRUE)))
      n <- length(x)
      pn <- private$pw2pn(m,parvect,stationary=stationary)
      foo <- pn$delta*dpois(x[1],pn$lambda)
      sumfoo <- sum(foo)
      lscale <- log(sumfoo)
      foo <- foo/sumfoo
      for (i in 2:n){
        if(!is.na(x[i])){
          P <- dpois(x[i],pn$lambda)
        }
        else {
          P <-rep(1,m)
        }
        foo <- foo %*% pn$gamma*P
        sumfoo <- sum(foo)
        lscale <- lscale+log(sumfoo)
        foo <- foo/sumfoo
      }
      mllk <- -lscale
      return(mllk)
    },
    # Computing the MLEs, given starting values for the natural parameters
    mle = function(delta = NULL){
      x <- private$x
      m <- private$m
      lambda0 <- private$lambda
      gamma0 <- private$gamma
      stationary <- private$stationary
      parvect0 <- private$pn2pw(m,lambda0,gamma0,delta,stationary=stationary)
      mod <- nlm(private$mllk,parvect0,x=x,m=m,stationary=stationary)
      pn <- private$pw2pn(m=m,mod$estimate,stationary=stationary)
      mllk <- mod$minimum
      np <- length(parvect0)
      AIC <- 2*(mllk+np)
      n <- sum(!is.na(x))
      BIC <- 2*mllk+np*log(n)
      return(list(m=m,lambda=pn$lambda,gamma=pn$gamma,delta=pn$delta,code=mod$code,
                  mllk=mllk,AIC=AIC,BIC=BIC))
    },
    # Forecast probabilities
    # Note that the output 'dxf' is a matrix.
    forecast_ = function(xf){
      h <- private$h
      x <- private$x
      mod <- private$fitted_model
      n <- length(x)
      nxf <- length(xf)
      dxf <- matrix(0,nrow=h,ncol=nxf)
      foo <- mod$delta*dpois(x[1],mod$lambda)
      sumfoo <- sum(foo)
      lscale <- log(sumfoo)
      foo <- foo/sumfoo
      for (i in 2:n){
        foo <- foo%*%mod$gamma*dpois(x[i],mod$lambda)
        sumfoo <- sum(foo)
        lscale <- lscale+log(sumfoo)
        foo <- foo/sumfoo
      }
      for (i in 1:h){
        foo <- foo%*%mod$gamma

        for (j in 1:mod$m){
          dxf[i,] <- dxf[i,] + foo[j]*dpois(xf,mod$lambda[j])
        }
      }

      return(dxf)
    }
  )
)
