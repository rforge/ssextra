\name{drawTreePop}
\alias{drawTreePop}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Generate a Synthetic Northern Hardoods Population

}

\description{

This routine will generate a simulated northern hardwoods population to
the user's specifications where diameters are drawn from a
three-parameter Weibull distribution. The heights and spatial process
are as described below. 

}

\usage{
drawTreePop(tract,
            solidTypes = c(1.5, 3),
            topDiams = c(0, 0),
            B = 80,
            hgt.sd = 6,
            a = 4,
            b = 8,
            c = 2,
            inhibitDist = 3,
            showPlot = TRUE,
            startSeed = 144,
            runQuiet = FALSE,
            ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{tract}{ A \code{"\linkS4class{Tract}"} or
           \code{"\linkS4class{bufferedTract}"} object in either
           \sQuote{English} or \sQuote{metric}. The units for this
           object specify the units for all other arguments, so be
           careful. }

  \item{solidTypes}{ See the argument of the same name in
	\code{\link[sampSurf:standingTree]{standingTree}}. }

  \item{topDiams}{ See the argument of the same name in
	\code{\link[sampSurf:standingTree]{standingTree}}. }

  \item{B}{ The basal area per acre, always in \sQuote{English}. }

  \item{hgt.sd}{ Standard deviation for height for \code{rnorm(0,
	  hgt.sd)} perturbations always in \sQuote{English}. }

  \item{a}{ Weibull location parameter always in \sQuote{English}. }

  \item{b}{ Weibull scale parameter always in \sQuote{English}. }

  \item{c}{ Weibull shape parameter always in \sQuote{English}. }

  \item{inhibitDist}{ This is the inhibition distance and it is in the
	same units as \code{tract}. Please see
	\code{\link[spatial:SSI]{SSI}} for more details. }

  \item{showPlot}{ \code{TRUE}: plot the sampled dbh distribution
     against the Weibull; FALSE: no plot. }

  \item{startSeed}{ A scalar specifying a random number seed
	    for the Weibull draws, spatial locations, &c. See
	    \code{\link[sampSurf:initRandomSeed]{initRandomSeed}} for
	    details. }

  \item{runQuiet}{ \code{TRUE}: no feedback; \code{FALSE}:
	    feedback if available. }

  \item{\dots}{ Further arguments passed on to
     \code{\link{sampleTrees}}. }

} %arguments

\details{

This routine will set up a reasonable simulated population of trees for
northern hardwoods based on the Weibull diameter distribution parameters
and height equation from Fast and Ducey (2011).

The spatial coordinates are drawn within the tract using an inhibition
process; the \pkg{spatial} package code is used here because it is faster
than the corresponding routine in the \pkg{spatstats} package.

The input units are always \sQuote{English} for everything but the
tract, which can be either; the output units are the same as that of the
tract. Therefore, if the tract is metric, then the returned tree list is
also.

}

\value{
 A data frame with the synthetic trees that can be cast to \code{\link{standingTrees}}.
}

\references{

Fast, A. J. and M. J. Ducey. 2011. Height-diameter equations for select New Hampshire tree species.
\emph{Northern Journal of Applied Forestry}, 28(3):157--160.

}

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

\note{

  Note especially the the basal area, height perturbations, and Weibull
  parameters are always entered in \sQuote{English} units, regardless of
  whether the \dQuote{Tract} object (and thus, the output tree data
  frame) is in metric. This is for people (like the author) that still
  have English as their primary units-of-measure.
  
} %note

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sampleTrees}}
}

\examples{
#
# a simple example...
#
\dontrun{
tr = Tract(c(x=100, y=100), cellSize=1)
btr = bufferedTract(10, tr)
s.trees = drawTreePop(btr, B=120) #BA = 27.5/ha
head(s.trees)
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
