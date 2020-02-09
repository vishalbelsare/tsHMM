\name{HMM.forecast}
\alias{HMM.forecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HMM.forecast(xf, h = 1, x, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xf}{
%%     ~~Describe \code{xf} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mod}{
%%     ~~Describe \code{mod} here~~
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
function (xf, h = 1, x, mod)
{
    n <- length(x)
    nxf <- length(xf)
    dxf <- matrix(0, nrow = h, ncol = nxf)
    foo <- mod$delta * dpois(x[1], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo/sumfoo
    for (i in 2:n) {
        foo <- foo \%*\% mod$gamma * dpois(x[i], mod$lambda)
        sumfoo <- sum(foo)
        lscale <- lscale + log(sumfoo)
        foo <- foo/sumfoo
    }
    for (i in 1:h) {
        foo <- foo \%*\% mod$gamma
        for (j in 1:mod$m) {
            dxf[i, ] <- dxf[i, ] + foo[j] * dpois(xf, mod$lambda[j])
        }
    }
    return(dxf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
