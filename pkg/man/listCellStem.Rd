\name{listCellStem}
\alias{listCellStem}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Enumerate Cell Numbers Within Inclusion Zones
  
}

\description{
  
  Perhaps a somewhat poor choice of name for this routine, it finds the
  cell numbers that are within each stem's inclusion zone, and vice versa
  (not within the actual stem!).

}

\usage{
listCellStem(ss, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ss}{ An object of class  \code{"\linkS4class{sampSurf}"}. }

  \item{\dots}{ Gobbled presently. }

} %arguments

\details{

  This routine forms the \dQuote{guts} of the
  \code{"\linkS4class{ssCellStemList}"} class constructor. For every
  inclusion zone, it determines which cell numbers fall within the zone
  (\code{cellsInStems}). Conversely, it will find what would be the
  \sQuote{inverse} of this as the stem numbers whose inclusion zones
  overlap each cell in the tract (\code{stemsInCells}). This routine is
  essential for determining the tree-based quantities like individual
  tree \acronym{VBAR}s that occur on a given sample point (grid
  cell). This allows us to aggregate tree-wise totals, &c., on a
  point-wise basis.

  This may be run \dQuote{stand-alone}, but it is probably better to
  call it as part of building an object of class
  \code{"\linkS4class{ssCellStemList}"}. In this way, any helper
  functions, &c., that are developed for that class in the future would
  be available for used on the results. In addition, that class contains
  more slots that may be helpful.
  
}

\value{
  
  A list with\ldots
  
\item{cellsInStems}{ A list by tree IDs (see
      \code{\link[sampSurf:getID]{getID}} for information on stem IDs)
      of the grid cell numbers within each tree's inclusion zone.  }

\item{stemsInCells}{ A list of all of the grid cell numbers with elements the IDs for 
      the stems whose inclusion zones overlap the individual 
      cells. Note that background cells are excluded since the list is  constructed from 
      \code{cellsInStems}. }

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

\seealso{
\code{"\linkS4class{ssCellStemList}"}, \code{"\linkS4class{ssBigBAF}"} %, \code{\link{wtVbars}}
}

\examples{
#
# if nhAll310.bb is an "ssBigBAF" object, then...
#
\dontrun{
#the following is for the big BAF surfaces...
csl.bb = listCellStem(nhAll310.bb@ss.bb.vol)
#and for the count surfaces...
csl.ct = listCellStem(nhAll310.bb@ss.ct.ba)
} %dontrun
} %examples

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
