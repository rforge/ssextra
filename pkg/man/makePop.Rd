\name{makePop}
\alias{makePop}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

Create a Northern Hardwoods Tract and Tree Population for Simulation
  
}

\description{

  This routine calls \code{\link{initTract}} to create a
  \code{"\linkS4class{bufferedTract}"} object. It then creates a tree
  population using \code{\link{drawTreePop}} within the
  \dQuote{bufferedTract}. Finally, it uses that tree population to
  generate a \code{\link{standingTrees}} object.
  
}

\usage{
makePop(extents = c(x = 150, y = 150),
        cellSize = 1,
        bufferWidth = 14,
        units = "metric",
        description = "Tract and northern hardwoods population",
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

  \item{description}{ A character description of the
	\dQuote{bufferedTract} object. }

  \item{\dots}{ Further arguments passed on to
     \code{\link{drawTreePop}}. }

} %arguments

\details{

As noted above, this is a simple way to create a \dQuote{bufferedTract}
object, along with a \code{\link{standingTrees}} object that can be used
in \code{\link{sampSurf}} simulations.
  
}

\value{

  A list (invisibly) with\ldots\cr
  
\item{btr}{ The \dQuote{bufferedTract} object.}

\item{strees}{ The \code{\link{standingTrees}} object.}

\item{trees}{ The corresponding data frame for the tree population.}
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

%% }

\examples{
#
# create a tract with defaults and plot it...
#
\dontrun{
pop = makePop()
plot(pop$btr, gridCol='grey90')
plot(pop$strees,add=TRUE)
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
