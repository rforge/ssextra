\name{ratioVariance}
\alias{ratioVariance}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Volume to Basal Area Ratio Variance
}

\description{

  This routine will calculate the volume to basal area ratio variance
  estimate from an individual Monte Carlo simulation replicate and
  compare it with the value from the \code{\link{monte}} run. For more
  details see the User's Guide.
  
} %description

\usage{
ratioVariance(ssBB,
              ssBB.mc,
              n = 25,
              mcSample = 1,
              stype = "bigBAF",
              runQuiet = FALSE,
              ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ssBB}{ An object of class \dQuote{\code{\linkS4class{ssBigBAF}}}. }

  \item{ssBB.mc}{ An object of class \dQuote{\code{\linkS4class{monteBigBAF}}}. }

  \item{n}{ The sample size, \eqn{n}, desired from the simulations. }

  \item{mcSample}{ The Monte Carlo sample replicate number for the specified sample size,
      \eqn{n}. }

  \item{stype}{  \sQuote{bigBAF} for tree stats on the big \acronym{BAF} surfaces; \sQuote{count} for
             the count surfaces. }

  \item{runQuiet}{ \code{TRUE}: run without any feedback; \code{FALSE}:
	some feedback desired. }

  \item{\dots}{ Gobbled. }

} %arguments

\details{

  This routine will calculate the ratio variance for the volume-basal area ratio in (8.34)
  of Gregoire and Valentine (2008). The variance is calculated via their equations (8.42)
  & (8.43); both being found on p. 258--259. The estimate is calculated for a single Monte
  Carlo sample at a given sample size.

} %details

\value{

  \item{msBB.mc}{The results of a call to \code{\link{monteStatsBB}} for \code{n} and
      \code{mcSample}. }

  \item{mtsBB.mc}{ The results of a call to \code{\link{monteTreeStatsBB}} for \code{n}
      and \code{mcSample}.}

  \item{n}{ The requested sample size. }

  \item{mcSample}{ The requested Monte Carlo sample replicate number. }

  \item{df}{ A data frame with the variance of the mean and standard errors for the ratio
      and normal theory estimators. }

  \item{means}{ A vector of means. }
  
}

\references{

  See the vignette for this package---\dQuote{A User's Guide to the ssExtra Package}---for
  an example and more information on the variance estimator and its use in the Monte Carlo
  simulations.
  
}

\author{
Jeffrey H. Gove
}

\note{

  Note that all of what is in this routine is also in the
  \code{"\linkS4class{monteBigBAF}"} constructor, but not through the use of this routine,
  which is simply made available as an easy way to calculate the variances for individual
  \acronym{MC} draws if desired, outside the framework of the constructor. The code in
  this routine is probably easier to understand than looking at the code in the
  constructor---specifically the requisite variances are calculated in
  \code{\link{monteStatsBB}}. It also demonstrates how one can extract the ratio variance
  estimate as well as a number of other quantities from the simulation results. This was
  originally the prototype code for the ratio variance calculations before it was added to
  the simulations through \code{monteStatsBB}.
  
} %note


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
   \code{\link{monteStatsBB}}
 }
 
\examples{
# Please see the user's guide for an example as it requires the output
# from a monte simulation run and is too lengthy to include here.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
