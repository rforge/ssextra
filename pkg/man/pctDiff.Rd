\name{pctDiff}
\alias{pctDiff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Utility to Calculate Percent Difference in Two Quantities
}

\description{

  This is a very simple little function to calculate the percentage
  difference between an \sQuote{estimate} and the \sQuote{true} value.
  
}

\usage{
pctDiff(est,
        true,
        ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{est}{ The \sQuote{estimate} value, which can be a vector, to test
	against the \code{true} value(s) below. }

  \item{true}{ The \sQuote{true} value, which can be a vector, to test
      against the \code{estimate} value(s) above. }
  
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}

} %arguments

\details{

  The simple percent difference: \eqn{(est-true)/true*100}.
  
} 

\value{

  The result of the above computation, which can be a vector or scalar.
  
} %value

%% \references{
%% %% ~put references to the literature/web site here ~
%% }

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }

\examples{
pctDiff(8, 10)
pctDiff(c(8, 12), 10)
pctDiff(c(8, 12), c(10, 6))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
