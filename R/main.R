#'@description Transforming natural parameters to working
#' @export
HMM.pn2pw <- function(m, lambda, gamma, delta = NULL, stationary = TRUE){

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

}



#'@description Transforming working parameters to natural
#' @export
HMM.pw2pn <- function(m,parvect,stationary=TRUE){

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

}



#'@description Computing minus the log-likelihood from the working parameters
#' @export
HMM.mllk <- function(parvect,x,m,stationary=TRUE,...){

  if(m==1) return(-sum(dpois(x,exp(parvect),log=TRUE)))
  n <- length(x)

  pn <- HMM.pw2pn(m,parvect,stationary=stationary)

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
}



#'@description Computing the MLEs, given starting values for the natural parameters
#' @export
HMM.mle <- function(x,m,lambda0,gamma0,delta=NULL,stationary=TRUE,...){

  parvect0 <- HMM.pn2pw(m,lambda0,gamma0,delta,stationary=stationary)

  mod <- nlm(HMM.mllk,parvect0,x=x,m=m,stationary=stationary)

  pn <- HMM.pw2pn(m=m,mod$estimate,stationary=stationary)

  mllk <- mod$minimum

  np <- length(parvect0)

  AIC <- 2*(mllk+np)

  n <- sum(!is.na(x))

  BIC <- 2*mllk+np*log(n)

  return(list(m=m,lambda=pn$lambda,gamma=pn$gamma,delta=pn$delta,code=mod$code,
              mllk=mllk,AIC=AIC,BIC=BIC))
}




#'@description Generating a sample
#' @export
HMM.generate_sample  <- function(ns,mod){
  mvect <- 1:mod$m

  state <- numeric(ns)

  state[1] <- sample(mvect,1,prob=mod$delta)

  for (i in 2:ns){
    state[i] <- sample(mvect,1,prob=mod$gamma[state[i-1],])
  }

  x <- rpois(ns,lambda=mod$lambda[state])

  return(x)
}

#'@description Global decoding by the Viterbi algorithm
#' @export
HMM.viterbi <-function(x,mod){

  n <- length(x)

  xi <- matrix(0,n,mod$m)

  foo <- mod$delta*dpois(x[1],mod$lambda)

  xi[1,] <- foo/sum(foo)
  for (i in 2:n){

    foo<-apply(xi[i-1,]*mod$gamma,2,max)*dpois(x[i],mod$lambda)

    xi[i,] <- foo/sum(foo)
  }
  iv<-numeric(n)

  iv[n] <-which.max(xi[n,])

  for (i in (n-1):1){
    iv[i] <- which.max(mod$gamma[,iv[i+1]]*xi[i,])
  }
  return(iv)
}




#'@description Computing log(forward probabilities)
#' @export
HMM.lforward <- function(x,mod){

  n <- length(x)
  lalpha <- matrix(NA,mod$m,n)
  foo <- mod$delta*dpois(x[1],mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- lscale+log(foo)

  for (i in 2:n){
    foo <- foo%*%mod$gamma*dpois(x[i],mod$lambda)
    sumfoo  <- sum(foo)
    lscale  <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo)+lscale
  }

  return(lalpha)
}



#'@description Computing log(backward probabilities)
#' @export
HMM.lbackward <- function(x,mod){
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA,m,n)
  lbeta[,n] <- rep(0,m)
  foo <- rep(1/m,m)
  lscale <- log(m)
  for (i in (n-1):1){
    foo <- mod$gamma%*%(dpois(x[i+1],mod$lambda)*foo)
    lbeta[,i] <- log(foo)+lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale+log(sumfoo)
  }
  return(lbeta)
}




#'@description Conditional probabilities
#==Conditional probability that observation at time t equals
#  xc, given all observations other than that at time t.
#  Note: xc is a vector and the result (dxc) is a matrix.
#' @export
HMM.conditional <- function(xc,x,mod){
  n <- length(x)

  m <- mod$m

  nxc <- length(xc)
  dxc <- matrix(NA,nrow=nxc,ncol=n)
  Px  <- matrix(NA,nrow=m,ncol=nxc)

  for (j in 1:nxc) {
    Px[,j] <-dpois(xc[j],mod$lambda)
  }

  la <- HMM.lforward(x,mod)
  lb <- HMM.lbackward(x,mod)
  la <- cbind(log(mod$delta),la)
  lafact <- apply(la,2,max)
  lbfact <- apply(lb,2,max)

  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*%mod$gamma)*exp(lb[,i]-lbfact[i])
    foo <- foo/sum(foo)
    dxc[,i] <- foo%*%Px
  }
  return(dxc)
}



#'@description Pseudo-residuals
#' @export
HMM.pseudo_residuals <- function(x,mod){
  n <- length(x)
  cdists <- HMM.conditional(xc=0:max(x),x,mod)
  cumdists <- rbind(rep(0,n),apply(cdists,2,cumsum))
  ulo <- uhi <- rep(NA,n)

  for (i in 1:n) {
    ulo[i]  <- cumdists[x[i]+1,i]
    uhi[i]  <- cumdists[x[i]+2,i]
  }
  umi <- 0.5*(ulo+uhi)
  npsr <- qnorm(rbind(ulo,umi,uhi))
  return(npsr)
}




#'@description State probabilities
#' @export
HMM.state_probs <- function(x,mod){
  n <- length(x)
  la <- HMM.lforward(x,mod)
  lb <- HMM.lbackward(x,mod)
  c <- max(la[,n])
  llk <- c+log(sum(exp(la[,n]-c)))
  stateprobs <- matrix(NA,ncol=n,nrow=mod$m)
  for (i in 1:n) {
    stateprobs[,i]<-exp(la[,i]+lb[,i]-llk)
  }

  return(stateprobs)
}


#'@description State prediction
# Note that state output 'statepreds' is a matrix even if h=1.
#' @export
HMM.state_prediction <- function(h=1,x,mod){
  n <- length(x)
  la <- HMM.lforward(x,mod)
  c <- max(la[,n])
  llk <- c+log(sum(exp(la[,n]-c)))
  statepreds <- matrix(NA,ncol=h,nrow=mod$m)
  foo <- exp(la[,n]-llk)

  for (i in 1:h){
    foo<-foo%*%mod$gamma
    statepreds[,i]<-foo
  }

  return(statepreds)
}


#'@description Local decoding
#' @export
HMM.local_decoding <- function(x,mod){
  n <- length(x)
  stateprobs <- HMM.state_probs(x,mod)
  ild <- rep(NA,n)
  for (i in 1:n) {
    ild[i]<-which.max(stateprobs[,i])
  }

  return(ild)
}



#'@description Forecast probabilities
# Note that the output 'dxf' is a matrix.
#' @export
HMM.forecast <- function(xf,h=1,x,mod){
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
