\name{ALFDA}
\alias{ALFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Affinity based Local Fisher Discriminant Analysis with network weights
%%  ~~function to do ... ~~
}
\description{Affinity based Local Fisher Discriminant Analysis with network weights
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ALFDA(genmat, y, r, kaf = 10, sigma = 0.5, 
metric = c("weighted", "orthonormalized", "plain"), 
knn = 6, reg = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genmat}{Genotype matrix or dataframe with markers p in columns and individuals n in rows
%%     ~~Describe \code{genmat} here~~
}
  \item{y}{n dimensional vector of class labels
%%     ~~Describe \code{y} here~~
}
  \item{r}{dimensionality of reduced space 
%%     ~~Describe \code{r} here~~
}
  \item{kaf}{ The number of nearest neighbours used in affinity matrix computation
%%     ~~Describe \code{kaf} here~~
}
  \item{sigma}{scaling parameter for affinity matrix
%%     ~~Describe \code{sigma} here~~
}
  \item{metric}{type of metric in the embedding space (default: 'weighted') 'weighted' - weighted eigenvectors 'orthonormalized' - orthonormalized 'plain' - raw eigenvectors
%%     ~~Describe \code{metric} here~~
}
  \item{knn}{parameter used in local scaling method (default: 6)
%%     ~~Describe \code{knn} here~~
}
  \item{reg}{regularization parameter (default: 0.001)
%%     ~~Describe \code{reg} here~~
}
}
\details{ Carrying out Affinity based Local Fisher Discriminant Analysis on genetic data and beyond
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Reduced features and genetic distance 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Sugiyama, M (2007). contain implementation Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}
\keyword{ALFDA}% use one of  RShowDoc("KEYWORDS")

