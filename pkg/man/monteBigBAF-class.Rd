\name{monteBigBAF-class}
\Rdversion{1.1}
\docType{class}
\alias{monteBigBAF-class}
%%\alias{hist,monteBigBAF-method}
\alias{show,monteBigBAF-method}
\alias{summary,monteBigBAF-method}

\title{Class \code{"monteBigBAF"}}

\description{
  
  This class defines a subclass of the virtual \code{"\linkS4class{monteDoubleSampling}"}
  for use with Monte Carlo simulations of Big \acronym{BAF} sampling.
  
}

\section{Objects from the Class}{
  
Objects can be created by calls to the constructor: \code{\link{monte}}
by passing a signature object of class
\code{"\linkS4class{ssBigBAF}"}. Using \code{new("monteBigBAF", ...)}
directly is \emph{discouraged} because of the class' complexity.

}

\section{Slots}{

  In addition to those slots that are defined within the
  \code{"\linkS4class{monteDoubleSampling}"}, the following are
  available\ldots

  \describe{
	
    \item{\code{mcSamples}:}{Object of class \code{"numeric"}: The
    number of Monte Carlo replications (samples). }
	
    \item{\code{n}:}{Object of class \code{"numeric"}: A vector of
    samples sizes to be drawn from the population (surfaces) and
    replicated \code{mcSamples} times. }
	
    \item{\code{fpc}:}{Object of class \code{"numeric"}: Finite
    population correction factor. }
	
    \item{\code{alpha}:}{Object of class \code{"numeric"}: The
    two-tailed \eqn{alpha} level for normal theory confidence interval
    construction. }
	
    \item{\code{replace}:}{Object of class \code{"logical"}:
    \code{TRUE}: each Monte Carlo replicate of size \code{n} is drawn
    with replacement; \code{FALSE}: samples are drawn withoug
    replacement from the surfaces.}
	
    \item{\code{ranSeed}:}{Object of class \code{"numeric"}: The
    starting seed for the random number generator. }
	
    \item{\code{t.values}:}{Object of class \code{"numeric"}: Student's
      \emph{t} values for each sample size \code{n} with two-tailed
      \eqn{alpha}-level \code{alpha}.}
	
    \item{\code{boot}:}{Object of class \code{"logical"}: \code{TRUE}:
    include jackknife and bootstrap estimates; \code{FALSE}: do not
    include these. }
	
    \item{\code{numBSS}:}{Object of class \code{"numeric"}: The number
    of bootstrap samples to be drawn from each Monte Carlo replicate and
    sample size if \code{boot=TRUE}.}
	
    \item{\code{means}:}{Object of class \code{"list"}: A
     list of\code{mcSamples x length(n)} data frames containing
       means for volume and basal area, tree VBARs, &c. }
	 
    \item{\code{vars}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing the
    sample variances. }
	
    \item{\code{stDevs}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing the
    sample standard deviations. }
	
    \item{\code{varMeans}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)}  data frames containing the
    sample variances of the means; e.g., Delta Method, Goodman, &c. }
	
    \item{\code{stErrs}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing the
    sample standard errors of the means. }
	
    \item{\code{lowerCIs}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing the
    lower confidence points on volume. }
	
    \item{\code{upperCIs}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing the
    upper confidence points on volume. }
	
    \item{\code{caught}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing the
    percent catch statistics for the above confidence intervals. }
	
    \item{\code{otherVarms}:}{Object of class \code{"list"}: A
    list of \code{mcSamples x length(n)} data frames containing other
    variances of the means of interest. }
	
    \item{\code{n.tvbar}:}{Object of class \code{"list"}: A list of
    \code{mcSamples x length(n)} data frames containing the number of tree VBARS
    on each replicate for the count and big \acronym{BAF} samples.}
  
    \item{\code{corrs}:}{Object of class \code{"list"}: A list of
    \code{mcSamples x length(n)} data frames containing the
    \emph{approximate aggregate} correlations of various types.}
	
    \item{\code{covs}:}{Object of class \code{"list"}: A list of \code{mcSamples x length(n)}
    data frames containing the \emph{point-based delta method} paired covariances. There
    are no covariances comparable to those in the \code{corrs} slot stored since their
    usefulness is questionable.}
	
    \item{\code{gm.all}:}{Object of class \code{"list"}: The grand
    summary means of the above: \code{means, vars, stDevs, varMeans,
    stErrs, lowerCIS, upperCIs, caught, otherVarms, corrs}. }
	
    \item{\code{sm.all}:}{Object of class \code{"list"}: The sampling
    variances and standard errors of the \code{means}. }
	
    \item{\code{mc.samples}:}{Object of class \code{"list"}: sample size
    list of data frames of size \code{length(n) x mcSamples} (for each
    \code{n}) holding the actual cell numbers drawn in each Monte Carlo
    sample replicate. }
	
  }
}

\section{Extends}{
Class \code{"\linkS4class{monteDoubleSampling}"}, directly.
}

\section{Methods}{
  \describe{
	
    \item{hist}{\code{signature(x = "monteBigBAF")}: Generate histograms }
	
    \item{show}{\code{signature(object = "monteBigBAF")}: Display a
    summary of the object. }
	
    \item{summary}{\code{signature(object = "monteBigBAF")}: Summary
    method for the object. }
	
	 }
} %section
   
%% \references{
%% %%  ~~put references to the literature/web site here~~
%% }

\author{
Jeffrey H. Gove
}

\note{ This describes a long complicated object structure. The constructor is also a long
  routine and takes time to run. It is best not to do bootstrapping if the number of
  \acronym{MC} samples is large as the bootstrap code is all in R and is very slow.  }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  The parent (super) class: \code{\linkS4class{monteDoubleSampling}}
}

\examples{
showClass("monteBigBAF")
}

\keyword{classes}
