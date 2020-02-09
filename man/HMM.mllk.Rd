\name{HMM.mllk}
\alias{HMM.mllk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing minus the log-likelihood from the working parameters
}
\description{
Computing minus the log-likelihood from the working parameters
}
\usage{
HMM.mllk(parvect, x, m, stationary = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parvect}{
%%     ~~Describe \code{parvect} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (parvect, x, m, stationary = TRUE, ...)
{
    if (m == 1)
        return(-sum(dpois(x, exp(parvect), log = TRUE)))
    n <- length(x)
    pn <- HMM.pw2pn(m, parvect, stationary = stationary)
    foo <- pn$delta * dpois(x[1], pn$lambda)
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo/sumfoo
    for (i in 2:n) {
        if (!is.na(x[i])) {
            P <- dpois(x[i], pn$lambda)
        }
        else {
            P <- rep(1, m)
        }
        foo <- foo \%*\% pn$gamma * P
        sumfoo <- sum(foo)
        lscale <- lscale + log(sumfoo)
        foo <- foo/sumfoo
    }
    mllk <- -lscale
    return(mllk)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
