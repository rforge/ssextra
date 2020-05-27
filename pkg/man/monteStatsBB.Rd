\name{monteStatsBB}
\alias{monteStatsBB}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Calculate Point-wise Statistics for Big \acronym{BAF} Sampling
}

\description{
  
This routine will calculate the simple point-wise statistics (mean, var,
var of the mean, etc.) on the individual sample draws from each of the
four sampling surfaces that can occur in big \acronym{BAF} sampling (see, e.g.,
"\code{\linkS4class{ssBigBAF}}").

} %description
 
\usage{
monteStatsBB(sdx = NA,
             ssBB,
             runQuiet = TRUE,
             ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{sdx}{ A vector of the cell numbers in the sample for a given
	Monte Carlo draw.}

  \item{ssBB}{ An object of class \dQuote{\code{\linkS4class{ssBigBAF}}}.}

  \item{runQuiet}{ \code{TRUE}: run without any feedback; \code{FALSE}:
	some feedback desired. }

  \item{\dots}{ Just gobbled for now. }
  
} %arguments
%} %usage

\details{
  
This routine is meant to be called from within \code{\link{monte}} in
the Monte Carlo sampling for construction of
\dQuote{\code{\linkS4class{monteBigBAF}}} objects. However, given a
sample of points (actually raster cells) from the sampling surfaces, one
can call this routine directly to calculate sampling statistics over all
points in the sample. For example, each
\dQuote{\code{\linkS4class{monteBigBAF}}} object contains all of the sample
replicate draws in the \code{mc.samples} slot, so one could pass the
\eqn{i}{i}th such draw to this routine as desired.

}

\value{
  
  A list with\ldots
  
\item{samples}{ A data frame with the point/cell-based sample (total) values for
  the main attribute from each sampling surface.}

\item{stats.ct}{ A data frame of statistics for the count surfaces.}

\item{stats.bb}{ A data frame of statistics for the big \acronym{BAF} surfaces.}

} %value

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
\code{\link{monteTreeStatsBB}}
}

\examples{
#
# if nhAll310.mc is an object of class "monteBigBAF"
# and nhAll310.bb is the corresponding "ssBigBAF" object
# then the following retrieves the first sample of
# size n=10, and the associated sample statistics...
#
\dontrun{
ns = nhAll310.mc@mc.samples$n.10$mc.1
stats = monteStatsBB(ns, nhAll310.bb)
}
} %examples
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
