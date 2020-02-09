\name{HMM.conditional}
\alias{HMM.conditional}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional probabilities.
Conditional probability that observation at time t equals xc, given all observations other than that at time t. Note: xc is a vector and the result (dxc) is a matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HMM.conditional(xc, x, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xc}{
%%     ~~Describe \code{xc} here~~
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
Zucchini, W., MacDonald, I., Langrock, R. (2016). Hidden Markov Models for Time Series. New York: Chapman and Hall/CRC, https://doi.org/10.1201/b20790}
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
function (xc, x, mod)
{
    n <- length(x)
    m <- mod$m
    nxc <- length(xc)
    dxc <- matrix(NA, nrow = nxc, ncol = n)
    Px <- matrix(NA, nrow = m, ncol = nxc)
    for (j in 1:nxc) {
        Px[, j] <- dpois(xc[j], mod$lambda)
    }
    la <- HMM.lforward(x, mod)
    lb <- HMM.lbackward(x, mod)
    la <- cbind(log(mod$delta), la)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    for (i in 1:n) {
        foo <- (exp(la[, i] - lafact[i]) \%*\% mod$gamma) * exp(lb[,
            i] - lbfact[i])
        foo <- foo/sum(foo)
        dxc[, i] <- foo \%*\% Px
    }
    return(dxc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
