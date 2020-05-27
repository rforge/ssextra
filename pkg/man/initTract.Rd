\name{initTract}
\alias{initTract}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Initialize a \code{"\linkS4class{bufferedTract}"} }

\description{

  This routine will set up a \sQuote{skeleton}
  \code{"\linkS4class{bufferedTract}"} object for simulation.

}

\usage{
initTract(extents = c(x = 150, y = 150),
          cellSize = 1,
          bufferWidth = 14,
          units = "metric",
          ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{extents}{ The x,y extents relative from (0,0) for the tract. Be
     sure to include the names (x & y as above) in the specification of
     this vector. }

  \item{cellSize}{ The cell size in feet or meters. }

  \item{bufferWidth}{ The buffer width, which should be large enough to
     include the half-width of the largest inclusion zone for any
     sampling method tested. }

  \item{units}{ \sQuote{English} or \sQuote{metric} (default). }

  \item{\dots}{ Other arguments that might be passed along to the
           \code{"\linkS4class{Tract}"} or
           \code{"\linkS4class{bufferedTract}"} constructors. }

} %arguments

\details{

  This is a very simple method whose only purpose is to make it easy to
  generate an \dQuote{bufferedTract} object.

} %details

\value{
An object of class  \dQuote{bufferedTract}.
} %value

%% \references{
%% %% ~put references to the literature/web site here ~
%% }

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

\note{
Note that \code{\link{makePop}} calls this directly, so the arguments
  above can be passed through that routine if creating a synthetic
  northern hardwoods population (through a call to \code{\link{drawTreePop}}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Tract}}, \code{\link{bufferedTract}}, \code{\link{makePop}}, \code{\link{drawTreePop}}
}

\examples{
#
# the default is metric; create a tract in acres with default
# one foot resolution...
#
\dontrun{
bt.eng = initTract(bufferWidth=20, units='English')
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
