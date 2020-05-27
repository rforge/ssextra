\name{monte-methods}
\docType{methods}
\alias{monte-methods}
%% \alias{monte,montePop-method}
%% \alias{monte,numeric-method}
%% \alias{monte,sampSurf-method}
\alias{monte,monteBigBAF-method}
\alias{monte,ssBigBAF-method}

\title{ Methods for "\code{\linkS4class{monteDoubleSampling}}" Object Construction in Package \pkg{ssExtra} }

\description{

  The following methods for generic function \code{\link{monte}} will construct valid
  objects of class "\code{\linkS4class{monteDoubleSampling}}".  Please see
  \code{\link[sampSurf:monte-methods]{monte-methods}} in \pkg{sampSurf} for more details
  on other more general methods for the \code{monte} generic.

}

\section{Methods}{

  %% \emph{Please note} that the signature argument in each case specifies
  %% the \emph{population} from which samples should be drawn.
  
\describe{

% \item{\code{signature(object = "montePop")}}{ Please see: \code{\link[sampSurf:monte-methods]{monte-methods}}
% in \pkg{sampSurf} for more details on this method.
% }

% \item{\code{signature(object = "numeric")}}{ Please see: \code{\link[sampSurf:monte-methods]{monte-methods}}
% in \pkg{sampSurf} for more details on this method.
% }

% \item{\code{signature(object = "sampSurf")}}{ Please see: \code{\link[sampSurf:monte-methods]{monte-methods}}
% in \pkg{sampSurf} for more details on this method.
% }



%--------------------------- ssBigBAF ----------------------------------------------------
\item{\strong{\code{signature(object = "ssBigBAF")}}}{

   This constructor method is only for use with objects of class
   "\code{\linkS4class{ssBigBAF}}".  } %item
  
  \item{}{\bold{usage\ldots}
    \preformatted{monte(object,
      n = c(10),
      mcSamples = 1,                   #number of MC samples
      alpha = 0.05,                    #two-tailed alpha level
      replace = TRUE,                  #MC sample with[out] replacement
      boot = FALSE,
      numBSS = 100,                    #number of bootstrap samples
      description = 'Big BAF Sampling Monte Carlo Simulation',
      startSeed = 123,
      runQuiet = FALSE,
      sshh = FALSE,
      ... )}
	   
   \itemize{
		 
      \item{\code{object}:}{ An object of class "\code{\linkS4class{ssBigBAF}}". }

	  \item{\code{n}:}{ A vector of sample sizes, \eqn{n}{n}.}

	  \item{\code{mcSamples}:}{ The number of Monte Carlo sample
	  replicates to draw for each sample size.}

	  \item{\code{alpha}:}{ The two-tailed alpha level for confidence
        interval calculation. I.e., for the 95\% confidence level,
        specify \code{alpha=0.05}.}

	  \item{\code{replace}:}{ See \code{\link{sample}} --- for sampling
	  with or without replacement.}
	   
	  \item{\code{boot}:}{ \code{TRUE}: include bootstrap & jacknife intervals;
	  \code{FALSE}: do not resample.}
	   
	  \item{\code{numBSS}:}{ The number of bootstrap samples if
		\code{boot}.}
	   
	  \item{\code{description}:}{ A character description for the return
	  object.}
		 
	  \item{\code{startSeed}:}{ A scalar specifying a random number seed
	    for the Monte Carlo draws. See
	    \code{\link[sampSurf:initRandomSeed]{initRandomSeed}} for details. }
		 
	  \item{\code{runQuiet}:}{ \code{TRUE}: no feedback; \code{FALSE}:
	    feedback if available. }
		 
	  \item{\code{sshh}:}{ \code{TRUE}: no extra feedback; \code{FALSE}:
	    allows a little feedback, even if \code{runQuiet = TRUE}. }
		
	  \item{... :}{ Other arguments to be passed along perhaps
	  eventually, right now they are \dQuote{gobbled.} }
	  
	} %itemize
 } %ssbigBAF item 
  
%} %ssBigBAF method






%--------------------------- monteBigBAF ------------------------------------------------

\item{\strong{\code{signature(object = "monteBigBAF")}}}{

    This constructor method is only for use with objects of class
    "\code{\linkS4class{monteBigBAF}}". It will accept an object from the above
    constructor and calculate the \emph{point-based delta method} (\acronym{PBDM})
    statistics for that object. The \acronym{PBDM} was not available at the time of
    writing of the above constructor---which is quite long---so this \dQuote{extention}
    constructor was written for this purpose.  } %item
  
  \item{}{\bold{usage\ldots}
    \preformatted{monte(object,
      ssBB = NA,
      runQuiet = FALSE,
      ... )}
	   
	\itemize{
		 
      \item{\code{object}:}{ An object of class "\code{\linkS4class{monteBigBAF}}". }

	  \item{\code{ssBB}:}{ An object of class "\code{\linkS4class{ssBigBAF}}". }
		 
	  \item{\code{runQuiet}:}{ \code{TRUE}: no feedback; \code{FALSE}:
	    feedback if available. }
		
	  \item{... :}{ Other arguments to be passed along perhaps
	  eventually, right now they are \dQuote{gobbled.} }
	  
	} %itemize
 } %montebigBAF item 
  
%} %monteBigBAF method






} %describe
} %section
\keyword{methods}
%\keyword{ ~~ other possible keyword(s) ~~ }
