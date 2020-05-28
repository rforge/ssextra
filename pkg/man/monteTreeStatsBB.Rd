\name{monteTreeStatsBB}
\alias{monteTreeStatsBB}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Calculate Tree-wise Statistics for Big \acronym{BAF} Sampling
}

\description{
  
This routine will calculate the simple \emph{tree}-wise statistics
(mean, var, var of the mean, correlation) on the individual sample draws from
each of the four sampling surfaces that can occur in big \acronym{BAF} sampling
(see, e.g., "\code{\linkS4class{ssBigBAF}}").

}

\usage{
monteTreeStatsBB(sdx,
                 ssBB,
                 stype = c("bigBAF", "count"),
                 ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{sdx}{ The sample cell numbers (the sample point index); note that
           it will be renamed to \code{cellNum} in the data frames
           exported here. }

  \item{ssBB}{ An object of class \dQuote{\code{\linkS4class{ssBigBAF}}}. }

  \item{stype}{ \sQuote{bigBAF} for tree stats on the big \acronym{BAF} surfaces; \sQuote{count} for
             the count surfaces. }

  \item{\dots}{  Just gobbled for now. }
}

\details{

This routine is meant to be called from within \code{\link{monte}} in
the Monte Carlo sampling for construction of
\dQuote{\code{\linkS4class{monteBigBAF}}} objects. However, given a
sample of points (actually raster cells) from the sampling surfaces, one
can call this routine directly to calculate sampling statistics over all
points in the sample. For example, each
\dQuote{\code{\linkS4class{monteBigBAF}}} object contains all of the
sample replicate draws in the \code{mc.samples} slot, so one could pass
the \eqn{i}{i}th such draw to this routine as desired.

It is important to note that the \emph{tree}-wise basal area results are
for the individual (and cell-wise aggregate) tree basal areas; they are
\emph{not} expanded. The tree (and cell-wise aggregate) \acronym{VBAR}s
are also returned. This is what distinguishes it from
\code{\link{monteStatsBB}}, which calculates summaries based on
\emph{point}-wise expanded \emph{totals}.

The routine is called once for the count surface summaries (\code{stype
  = 'count'}), and then again for the big \acronym{BAF} summaries
  (\code{'bigBAF'}) as desired.

}

\value{
  
  A list with\ldots
  
\item{df.trees}{ A tree-wise data frame with each row containing the
  results (tree basal area and \acronym{VBAR}) from an individual tree
  in the sample, sorted by point. It also includes zero-records for any
  background cells selected in the sample of points. }

\item{list.df}{ A list of data frames corresponding to \code{df.trees},
     where the individual data frames are \code{NULL} for packground
     points, or the trees measured in that cell/point; this list will
     have length \code{length(sdx)} (i.e., the sample size in terms of
     number of points). }
  
\item{df.cells}{ A cell-wise (point-wise) data frame where the trees
   have been aggregated (summed) by cell for tree-based (\emph{not}
   expanded) totals at each cell/point. }
  
\item{treeStats}{ A data frame with the summary stats over all trees in
  \code{df.trees}. }
  
\item{cellStats}{ A data frame with the summary stats over all
  cells/points in \code{df.cells}. }
  
\item{n.v}{ The number of trees in the sample (i.e., in
  \code{df.trees}). This is the number of sampled \acronym{VBAR}
  trees. }
  
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
\code{\link{monteStatsBB}}
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
stats = monteTreeStatsBB(ns, nhAll310.bb)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
