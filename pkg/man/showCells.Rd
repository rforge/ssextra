\name{showCells}
\alias{showCells}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  
  Highlight desired cells on a Sampling Surface
  
}

\description{

  This function allows one to display certain cells with a different
  background color and, optionally mark their cell centers. This is
  usedful when it might be desired to highlight certain cells on a
  sampling surface. Potential uses include showing where a subsample of
  cells are located on the surface. Also, by using
  \code{\link{listCellStem}}, one can offset individual overlap regions
  of inclusion zones.
  
}

\usage{
showCells(ss,
          cells,
          cell.col = "salmon2",
          alphaTrans = 0.2,
          showSS = TRUE,
          showPoints = TRUE,
          point.col = "gray25",
          ...)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{ss}{ An object of class \code{\linkS4class{sampSurf}}.  }
  
  \item{cells}{ A vector of cell numbers in the
    \code{\linkS4class{sampSurf}} object to be displayed differently as sample
    points. }

  \item{cell.col}{ The color for the cells to be displayed. }
  
  \item{alphaTrans}{ Transparencey: see
    \code{\link[sampSurf]{transparentColorBase}}. }

  \item{showSS}{ \code{TRUE}: first plot the sampSurf object;
     \code{FALSE}: assume we are adding to an existing figure. Note that
     the existing figure need \emph{not} be a sampling surface, but must
     be of the same extents as the sampling surface passed in the
     \code{ss} argument. }

  \item{showPoints}{ \code{TRUE}: show the sample points (cell centers);
	\code{FALSE}: display only the entire cell shading. }

  \item{point.col}{ The color for the sample points (cell centers). }

  \item{\dots}{ Passed on to \code{\link{image}} and \code{\link{points}}. }

} %\arguments

\details{

  This routine can be used to either plot a \dQuote{sampSurf} object and
  add a vector of sample points to it for display, or simply to add the
  sample points to an existing plot. Either way, the assumption is that
  we have sampled the full sampling surface in the sense of drawing a
  subsample of points from the full set of cells in the
  \code{\linkS4class{Tract}} slot of the object. For example, this would
  allow the display of sample points drawn randomly from a surface in
  one Monte Carlo replication from any of the methods associated with
  the \code{\link[sampSurf]{monte}}, including the method(s) contained
  within the current package (\pkg{ssExtra}).

  Note that with the default choices for \code{cell.col} and
  \code{alphaTrans} it can be very difficult to see the shading for the
  cells on larger tracts (rasters). In this case, decrease the
  transparency and choose a brighter color for the cells in order to
  display them more prominently. Similarly for \code{point.col}, though
  \code{alphaTrans} is not directly used for the points; however, one
  can pass \code{point.col} as a color that has been returned from
  \code{\link[sampSurf]{transparentColorBase}} to achieve the same
  effect.
  
}

\value{

  The function returns an object of class
  \code{\linkS4class{RasterLayer}} (invisibly) where \emph{only} those
  cells that were passed in \code{cells} are \emph{not} \code{NA}. The
  values for the cells are the cell numbers passed.
  
}

%% \references{
%% %% ~put references to the literature/web site here ~
%% }

\author{
Jeffrey H. Gove
}

%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[ssExtra:plot-methods]{plot}}
}

\examples{
\dontrun{
tr = Tract(c(x=50,y=50), cellSize=0.5)
btr = bufferedTract(10, tr)
ssSA = sampSurf(3, btr, iZone = 'sausageIZ', plotRadius=3,
       buttDiam=c(30,50), startSeed=1001)
plot(ssSA, axes=TRUE)
summary(ssSA)
#
#add the sample points...
#
showCells(ssSA, cells = sample(1:10000, 100), showSS=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
