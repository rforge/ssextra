\name{ssCellStemList-class}
\Rdversion{1.1}
\docType{class}
\alias{ssCellStemList-class}

\title{Class \code{"ssCellStemList"}}

\description{ This class is used in the simulation of big \acronym{BAF}
  sampling through the \code{"\linkS4class{ssBigBAF}"} class slots. It
  is a class that probably would not be used too much outside of the
  above or similar. It contains several lists and data frames that aid
  in working with tree-wise aggregates like tree \acronym{VBAR}s on
  individual cells/points.

} %description

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ssCellStemList",
  ...)}. However, it is simpler and less error-prone to use the constructor of the
same name: \code{\link{ssCellStemList}}.
}

\section{Slots}{
  
  \describe{
	
    \item{\code{description}:}{Object of class \code{"character"} --- A
     description of the object. }
	
    \item{\code{cellsInStems}:}{Object of class \code{"list"} --- A
      list of cell numbers w/in each individual stem's inclusion zone. }
	
    \item{\code{stemsInCells}:}{Object of class \code{"list"} --- A list
     of stem IDs whose inclusion zones overlap each individual
     cell---for all cells that contain inclusion zones. Note that
     background cells are excluded. }
	
    \item{\code{sumSIC}:}{Object of class \code{"data.frame"} --- A
    summary data frame based on \code{stemsInCells} with the number of
    cells containing \eqn{n}{n} stem inclusion zone overlaps. }
	
    \item{\code{sumCIS}:}{Object of class \code{"data.frame"} --- A
    summary data frame based on \code{cellsInStems} with the number of
    stem inclusion zones containing \eqn{n}{n} cells. This will normally
    be quite a bit longer than the \code{sumSIC} data frame. }
	
    \item{\code{cellIDs}:}{Object of class \code{"character"} --- A
     vector of the cell IDs in the \code{stemsInCells} list. }
	
    \item{\code{cellNums}:}{Object of class \code{"integer"} --- Integer
    version of \code{cellIDs}. }
	
    \item{\code{stemIDs}:}{Object of class \code{"character"} --- A
     vector of the stem IDs in the \code{cellsInStems} list. }
	
  } %describe
} %section

\section{Methods}{
There are currently no methods defined with class "ssCellStemList" in
the signature.
}

%% \references{
%% %%  ~~put references to the literature/web site here~~
%% }

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{listCellStem}}
}

\examples{
#
# if nhAll310.bb is an "ssBigBAF" object, then
#
\dontrun{
head(nhAll310.bb@csl.bb@cellsInStems)
length(nhAll310.bb@csl.bb@cellsInStems)
length(nhAll310.bb@csl.bb@stemsInCells)
} %dontrun
}
\keyword{classes}
