\name{createBBNH}
\alias{createBBNH}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Create a \dQuote{\code{\linkS4class{ssBigBAF}}} object for \dQuote{Big \acronym{BAF}}
  Sampling Northern Hardwood Simulations
  
}

\description{

  This routine is used for the specific purpose of creating a
  \dQuote{\code{\linkS4class{ssBigBAF}}} object for a synthetic northern hardwoods
  stand. The object created can then be used in Monte Carlo simulation experiments to look
  more closely at the mechanics of \dQuote{Big \acronym{BAF}} sampling.

}

\usage{
createBBNH(extents = c(x = 178, y = 178),
           cellSize = 1,
           bufferWidth = 18,
           units = "metric",
           baf.ct = 4,
           baf.bb = 10,
           startSeed = 355,
           ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{extents}{ The extents of raster cells in \code{x} and \code{y}
	(from an origin at \code{(0,0)}) for the four
	\dQuote{\code{\linkS4class{bufferedTract}}} objects that will be
	created to hold the different sampling surfaces. }

  \item{cellSize}{ The cell size in meters (for \code{units = "metric"})
  or feet (for \code{units = "English"}). }

\item{bufferWidth}{ The width of the buffer (see
    \dQuote{\code{\linkS4class{bufferedTract}}} for more information). Please note that
    this should be large enough to include the half-width of the largest inclusion zone
    for any sampling method tested; i.e., in the case of \dQuote{Big \acronym{BAF}}
    sampling, the largest inclusion zone for the \code{baf.ct} basal area factor. }

  \item{units}{ \dQuote{metric} or \dQuote{English}. }
  
  \item{baf.ct}{ The basal area factor for the count sample. }

  \item{baf.bb}{ The basal area factor for the volume sample. }

  \item{startSeed}{ A seed for the random number generator. See
	\code{\link{initRandomSeed}} for more details. }

  \item{\dots}{ Arguments that are passed on to \code{\link{makePop}},
	the \code{\link{sampSurf}} constructor and \code{\link{ssBigBAF}}
	constructor. Please see these routines for the different arguments
	that they accept. }

} %argument


\details{

  The default values in the arguments to this routine will duplicate the
  northern hardwood sampling surfaces for one pair of basal area factors
  that are used in Gove et al (2020). A quick look at the \R code in
  this routine will show the necessary steps to make
  \dQuote{\code{\linkS4class{ssBigBAF}}} objects for other forest
  types. There is also some discussion of this in the package user's
  guide vignette.
  
} %details

\value{

  An object of class \dQuote{ssBigBAF}, invisibly.

} %value

\references{

\describe{
%

\item{}{Gove, J. H., Gregoire,T. G., Ducey, M. J., and Lynch,
   T. B. 2020. A Note on the Estimation of Variance for Big BAF
   Sampling. \emph{Forest Ecosystems}, Submitted.}

}%describe

  
} %references

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

\note{
  
  It may take some quick calculations to determine the width of the
  buffer required to hold the largest inclusion zone. And it may have to
  be iterated as well. But it is quite simple really.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{initTract}},
  \code{\link{makePop}},
  \code{\link{drawTreePop}},
  and the \dQuote{\code{\linkS4class{ssBigBAF}}} class.
}

\examples{
#
# use the defaults...
#
\dontrun{
ssBB.nh = createBBNH()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
