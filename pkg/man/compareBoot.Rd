\name{compareBoot}
\alias{compareBoot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Compare Multiple Bootstrap Simulation Resuts from \pkg{bcaboot} and \pkg{boot} 
}

\description{

  This routine will run \code{\link{testBoot}} on each of \eqn{n+B}{n+B}
  combinations. Thus, it can be used to look at the differences between
  the two bca confidence interval routines,
  \pkg{bcaboot:}\code{\link{bcajack}} and
  \pkg{boot:}\code{\link{boot.ci}} for confidence interval calculation.

}

\usage{
compareBoot(n = c(100, 200),
            B = c(1000, 1200),
            tablePath = "",
            tableName = "",
            digits = c(0, 0, 0, rep(3, 6)),
            caption = "Bootstrap Confidence Interval Comparison",
            label = "tab:bootTable",
            ssshhh = FALSE,
            ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{n}{ A vector of sample sizes that are to be used for
  bootstrapping. Note that since this routine calls
  \code{\link{testBoot}}, the samples for each sample size are drawn
  from a standard normal distribution. }

  \item{B}{ A vector with the number of bootstrap sample replicates at
  each sample size. }

  \item{tablePath}{ The relative path to \code{tableName}; see below for
  detials. }

  \item{tableName}{ The file name for the \pkg{xtable} output table. Use
  \sQuote{} for no hardcopy table output; otherwise a \code{*.tex} file
  name. }

  \item{digits}{ The number of display digits for the
  \code{\link[xtable:xtable]{xtable}}. Please see the link for more details. }

  \item{caption}{ The caption for the \code{\link[xtable:xtable]{xtable}}. }

  \item{label}{ The \emph{LaTeX} label for the \code{\link[xtable:xtable]{xtable}}. }

  \item{ssshhh}{ \code{runQuiet} is set to \code{TRUE} in the call to
    \code{testBoot}, so this is a local version only for running
    quietly. }

  \item{\dots}{ Arguments passed to \code{testBoot}. }
  
} %arguments

\details{
  
This routine was written simply to facilitate comparison the two bca
confidence interval results in the different packages cited above. Note
that the results can be saved to a file by using the
\code{\link[xtable:xtable]{xtable}} function in package \pkg{xtable}.
  
} %details

\value{

  A list with\ldots
  
  \item{ df }{ A data frame with the results, which should be
	self-explanatory.}
  
  \item{ xtab }{ A \emph{LaTeX} table object of class \dQuote{xtable}. }

} %value

%% \references{
%% %% ~put references to the literature/web site here ~
%% } %references

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

\note{
  
  Neither this routine, nor \code{testBoot} are meant to be general
  purpose routines. They were written to explore the differences in
  calculated bca intervals for a specific population, that of a standard
  normal. One could certainly adapt these routines to other uses.

  Please see the note in \code{\link{testBoot}} concerning the potential
  error that can occur using \code{\link{boot.ci}}.
  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{testBoot}},
  \code{\link{bcaboot}},
  \code{\link{bcajack}},
  \code{\link{boot}},
  \code{\link{boot.ci}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
#
# simply...
#
\dontrun{
cmp.bs = compareBoot()
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
