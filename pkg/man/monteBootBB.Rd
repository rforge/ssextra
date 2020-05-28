\name{monteBootBB}
\alias{monteBootBB}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  
Calculate Bootstrap and Jackkife Statistics on an Individual Monte Carlo
Sample

}

\description{

This routine will perform bootstrap resampling as well as the jackknife
on a single Monte Carlo sample replicate drawn from the sampling
surfaces of "\code{\linkS4class{ssBigBAF}}" objects. This routine is
called from within the nested sample size--within--Monte Carlo loops in
the "\code{\linkS4class{monteBigBAF}}" constructor. However, it
certainly can be called in other situations where the correct
information is available. Among other things, it calculates the
\sQuote{\acronym{bca}} confidence intervals on the replicate sample.
  
}

\usage{
monteBootBB(samples,
            B = 100,
            alpha = 0.05,
            ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{samples}{ This must be a \code{samples} data frame calculated in
               a single Monte Carlo draw as returned from a call to
               \code{\link{monteStatsBB}} for a given sample size,
               \eqn{n}{n}. }

  \item{B}{ The number of bootstrap samples. }

  \item{alpha}{ The two-tailed \eqn{\alpha}{alpha} level, used to select
	the correct bootstrap \sQuote{\acronym{bca}} intervals for
	comparison. }

  \item{\dots}{ Gobbled for now. }

} %arguments

\details{

  This routine uses the \code{\link{bcajack}} routine from \pkg{bcaboot}
for construction of both bootstrap and jackknife confidence intervals
for a Monte Carlo sample from an "\code{\linkS4class{ssBigBAF}}"
object. The calculations in that package are all done in native \code{R}
code, and so can be slow when applied to many Monte Carlo replicates in
the context of big \acronym{BAF} sampling. Please see the \pkg{bcaboot}
package for more details, and also the reference below.

}

\value{

A list invisibly as described in the return from \code{\link{bcajack}}\ldots

\item{call}{ The call to \code{bcajack} as returned from
  \code{\link{match.call}} }

\item{lims}{ \dQuote{first column shows the estimated bca confidence
	limits\ldots} }

\item{stats}{ \dQuote{top line of stats shows 5 estimates\ldots} }

\item{B.mean}{ \dQuote{bootstrap sample size B, and the mean\ldots} }

\item{ustats}{ \dQuote{The bias-corrected estimator\ldots} }

\item{seed}{ \dQuote{The random number state for reproducibility} }

including the following as part of the list\ldots

\item{jackmean}{ The jackknife estimate of the mean. }

} %value

\references{

  Efron, B. and B. Narasimhan. 2018. The automatic construction of
  bootstrap confidence intervals. Technical Report No. 2018-07,
  Standford University Department of Statistics, PDF version: 
  \href{https://statweb.stanford.edu/~ckirby/brad/papers/}{here}.

}

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

%% \note{
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

Please also see the \code{\link{testBoot}} routine for some more
information on a comparison of the results from the \pkg{bcaboot}
package to those in the \pkg{boot} package.

}

\examples{
#
# If 'stats' is a return object from monteStatsBB, then
# the following will draw 200 bootstrap samples and also
# calculate the jackknife results...
#
\dontrun{
mcb = monteBootBB(stats$samples, 200)
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
