\name{HMM.generate_sample}
\alias{HMM.generate_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating a sample}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HMM.generate_sample(ns, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ns}{
%%     ~~Describe \code{ns} here~~
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
function (ns, mod)
{
    mvect <- 1:mod$m
    state <- numeric(ns)
    state[1] <- sample(mvect, 1, prob = mod$delta)
    for (i in 2:ns) {
        state[i] <- sample(mvect, 1, prob = mod$gamma[state[i -
            1], ])
    }
    x <- rpois(ns, lambda = mod$lambda[state])
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
