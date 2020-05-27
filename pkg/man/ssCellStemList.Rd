\name{ssCellStemList}
\alias{ssCellStemList}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Generate Objects of Class "\code{\linkS4class{ssCellStemList}}"
}

\description{

This generic function has one method. It is used as a constructor
function for objects that are of class
"\code{\linkS4class{ssCellStemList}}" and should be used in preference
to \code{\link{new}} to insure a valid object.

} %description

\usage{
ssCellStemList(ss, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{ss}{ An oject of class "\code{\linkS4class{sampSurf}}". }

  \item{\dots}{ Gobbled at present. }

} %arguments

\details{

Objects of class "\code{\linkS4class{ssCellStemList}}" are created using
\code{\link{ssCellStemList-methods}} by passing an argument of class
\dQuote{sampSurf}. The guts of object creation is handled by the
\code{\link{listCellStem}} function. A lot of work goes into the
creation of these objects and so it can take a little time on larger
\dQuote{sampSurf} objects (i.e., higher resolution and larger extents).

}

\value{
An object of class "\code{\linkS4class{ssCellStemList}}".
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
\code{\link{listCellStem}}
}

\examples{
#
# if ssBB is an object of class "ssBigBAF" then...
#
\dontrun{
csl.bb = ssCellStemList(ssBB@ss.bb.vol) #BB list
csl.ct = ssCellStemList(ssBB@ss.ct.ba)  #count list
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
