\name{HMM.state_prediction}
\alias{HMM.state_prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
State probabilities
}
\description{
State probabilities
}
\usage{
HMM.state_prediction(h = 1, x, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
Resul Akay
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
function (h = 1, x, mod)
{
    n <- length(x)
    la <- HMM.lforward(x, mod)
    c <- max(la[, n])
    llk <- c + log(sum(exp(la[, n] - c)))
    statepreds <- matrix(NA, ncol = h, nrow = mod$m)
    foo <- exp(la[, n] - llk)
    for (i in 1:h) {
        foo <- foo \%*\% mod$gamma
        statepreds[, i] <- foo
    }
    return(statepreds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
