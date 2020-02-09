\name{HMM.mle}
\alias{HMM.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HMM.mle(x, m, lambda0, gamma0, delta0 = NULL, stationary = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{lambda0}{
%%     ~~Describe \code{lambda0} here~~
}
  \item{gamma0}{
%%     ~~Describe \code{gamma0} here~~
}
  \item{delta0}{
%%     ~~Describe \code{delta0} here~~
}
  \item{stationary}{
%%     ~~Describe \code{stationary} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Zucchini, W., MacDonald, I., Langrock, R. (2016). Hidden Markov Models for Time Series. New York: Chapman and Hall/CRC, https://doi.org/10.1201/b20790
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, m, lambda0, gamma0, delta0 = NULL, stationary = TRUE,
    ...)
{
    parvect0 <- HMM.pn2pw(m, lambda0, gamma0, delta0, stationary = stationary)
    mod <- nlm(HMM.mllk, parvect0, x = x, m = m, stationary = stationary)
    pn <- HMM.pw2pn(m = m, mod$estimate, stationary = stationary)
    mllk <- mod$minimum
    np <- length(parvect0)
    AIC <- 2 * (mllk + np)
    n <- sum(!is.na(x))
    BIC <- 2 * mllk + np * log(n)
    return(list(m = m, lambda = pn$lambda, gamma = pn$gamma,
        delta = pn$delta, code = mod$code, mllk = mllk, AIC = AIC,
        BIC = BIC))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
