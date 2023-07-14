\name{Af}
\alias{Af}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the affinity matrix that represents the neighborhood graph of the individuals
%%  ~~function to do ... ~~
}
\description{
Compute the affinity matrix that represents the neighborhood graph of the individuals
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Af(pwD, K = 10, sigma = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pwD}{pairwise distance /similarity matrix 
%%     ~~Describe \code{pwD} here~~
}
  \item{K}{ The number of nearest neighbours
%%     ~~Describe \code{K} here~~
}
  \item{sigma}{ scaling parameter 
%%     ~~Describe \code{sigma} here~~
}
}
\details{ Compute the affinity matrix that represents the neighborhood graph of the individuals, the affinity matrix represents the weights between individuals
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A affinity matrix representing the weights between individuals 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
qinxinghu@gmail.com
%%  ~~who you are~~
}

\keyword{ Affinity matrix}% use one of  RShowDoc("KEYWORDS")

