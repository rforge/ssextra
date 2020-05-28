\name{testBoot}
\alias{testBoot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Test \sQuote{bca} Confidence Intervals from \pkg{bcaboot} and \pkg{boot}
}

\description{

  This routine will calculate bootstrap \sQuote{bca} confidence
  intervals from a sample of standard normal random variates. The output
  will facilitate comparison of the corresponding routines for
  \sQuote{bca} calculation from two popular packages: \pkg{bcaboot} and
  \pkg{boot}.

}

\usage{
testBoot(n = 1000,
         B = 1000,
         alpha = 0.05,
         mbm = FALSE,
         parallel = "no",
         ncpus = 3L,
         times = 10L,
         unit = "s",
         runQuiet = FALSE,
         startSeed = 245,
         ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{n}{ The number of standard normal random variates in the
	sample. }

  \item{B}{ The number of bootstrap sample replicates. }

  \item{alpha}{ The two-tailed \eqn{\alpha}{alpha} level for comparison
    with normal theory intervals. }

  \item{mbm}{ \code{TRUE}: run microbenchmark for \code{times};
	\code{FALSE}; no timing. }

  \item{parallel}{ See \code{\link{boot}}, it does not seem to work. }

  \item{ncpus}{ See \code{\link{boot}}. }

  \item{times}{ The number of \code{\link{microbenchmark}} replicates. }

  \item{unit}{ The units for the report from microbenchmark; \sQuote{s}
	= seconds. }

  \item{runQuiet}{ \code{TRUE}: sssshhh; \code{FALSE}: print results. }

  \item{startSeed}{ A scalar specifying a random number seed
	    for the Monte Carlo draws. See
	    \code{\link[sampSurf:initRandomSeed]{initRandomSeed}} for
	    details. }

  \item{\dots}{ Other arguments passed to \code{\link{boot}}. }

} %arguments

\details{

  Calculating the bootstrap \sQuote{bca} intervals adds a lot of extra
  time to the big \acronym{BAF} Monte Carlo simulations. I use the
  \pkg{boot} package in \code{\link[sampSurf:monte]{monte}}, and it
  seems to have worked well. But recently the \pkg{bcaboot} package
  appeared and I wanted to see if it had any advantage, and how close
  the calculated confidence intervals were between the two. The results
  are also compared against the normal theory intervals. I assumed that
  \code{\link{bcajack}} from \pkg{bcaboot} might be somewhat better
  since it uses the jackknife standard error (though I am not sure about
  \code{boot}).

  I was somewhat disappointed in the results, which differ, and seem to
  take a large number of bootstrap replicates to converge to the normal
  theory intervals (this actually never happened, they just tend that
  way). Also, \code{boot} can throw an error (see \sQuote{Note} section
  below) that I had not encountered before. Of course, the
  \code{\link{boot}} function has a lot of options, which I did not try,
  but can be passed using the \dots argument above if desire.

  In the end, I decided to use \code{\link{bcajack}} for the intervals
  in the version of \code{\link{monte}} for big \acronym{BAF} sampling
  in the \pkg{ssExtra} package.
  
} %details

\value{

A list invisibly with\ldots

\item{samp}{ The standard normal sample drawn.}

\item{res}{ The results from \code{\link{bcajack}}.}

\item{boot.samp}{ The results from \code{\link{boot}}.}

\item{boot.cis}{ The results from \code{\link{boot.ci}}.}

\item{mnm.en}{ The microbenchmark results from \code{bcajack}.}

\item{mbm.cr}{ The microbenchmark results from \code{boot} &
  \code{boot.ci}.}

\item{df.samp}{ A data frame with the normal theory results from the
  sample.}

\item{df}{ A data frame with the bootstrap results.}

\item{call}{ The results from \code{\link{match.call}}. }

} %value

%% \references{
%% %% ~put references to the literature/web site here ~
%% }

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

\note{

   Be careful with \code{\link{boot.ci}} as it requires \eqn{B>n}{B>n} evidently or you will get
   an error. One time I got an error because \code{B} was not a multiple of \code{n}.
   It seems that making \code{B = 2n} is the usual trick; e.g., see 
   \href{https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html}{here}.

   The error produced is\ldots
   \code{bca.ci(boot.out, conf, index[1L], L = L, t = t.o, t0 = t0.o,  : 
  estimated adjustment 'a' is NA}

} %note

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{compareBoot}},
  \code{\link{bcaboot}},
  \code{\link{bcajack}},
  \code{\link{boot}},
  \code{\link{boot.ci}}
}

\examples{
#
# simply...
#
\dontrun{
bs = testBoot()
} %dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
