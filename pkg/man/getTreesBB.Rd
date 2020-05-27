\name{getTreesBB}
\alias{getTreesBB}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Extract the Tree Population from a \code{"\linkS4class{sampSurf}"} object
}

\description{

This is a utility routine that will extract the tree population from an
object of class \dQuote{sampSurf}. It will return the population in two
different formats as noted below.
  
}

\usage{
getTreesBB(ssObj, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ssObj}{ An object of class \code{"\linkS4class{sampSurf}"}. }

  \item{\dots}{ Gobbled currently. }

} %arguments

\details{

This routine will return all of the necessary individual tree attributes
for recreating the population as a \code{"\linkS4class{standingTrees}"}
object, including the spatial locations for the trees. Of course, that
is not necessary, since the \dQuote{standingTrees} object is included in
the return. However, in addition to the usual information for such an
object, the data frame returned also includes individual tree
\acronym{VBAR}s, which makes it useful for big \acronym{BAF} sampling.

Recall that the tree population is independent of the estimate attribute
for the sampling surface. Therefore, any \dQuote{sampSurf} object can be
passed to extract the population of trees (see, e.g., the example below).
  
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

A list invisibly with\ldots

\item{trees}{ A data frame containing the full tree population. }

\item{strees}{ An object of class \code{"\linkS4class{standingTrees}"}
  corresponding to the \code{trees} data frame.}

\item{nTrees}{ The number of trees in both \code{trees} and
  \code{strees}.}

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
\code{\link{standingTrees-methods}}, \code{\link{sampleTrees}}
}

\examples{
#
# if ssBB is a "ssBigBAF" object, then, simply...
#
\dontrun{
tp.vol = getTreesBB(ssBB@ss.bb.vol) #from volume
tp.ba = getTreesBB(ssBB@ss.bb.ba)   #from basal area
identical(tp.vol, tp.ba)  #identical
} %dontrun
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
