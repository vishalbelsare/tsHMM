\name{HMM.viterbi}
\alias{HMM.viterbi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HMM.viterbi(x, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (x, mod)
{
    n <- length(x)
    xi <- matrix(0, n, mod$m)
    foo <- mod$delta * dpois(x[1], mod$lambda)
    xi[1, ] <- foo/sum(foo)
    for (i in 2:n) {
        foo <- apply(xi[i - 1, ] * mod$gamma, 2, max) * dpois(x[i],
            mod$lambda)
        xi[i, ] <- foo/sum(foo)
    }
    iv <- numeric(n)
    iv[n] <- which.max(xi[n, ])
    for (i in (n - 1):1) {
        iv[i] <- which.max(mod$gamma[, iv[i + 1]] * xi[i, ])
    }
    return(iv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
