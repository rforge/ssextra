\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ssBigBAF,missing-method}

\title{ Methods for Function \code{plot} in Package \pkg{ssExtra} }

\description{

These are simple methods to plot objects found in this package; e.g.,
\dQuote{ssBigBAF} objects, which have slots that are sampling
surfaces. Many more methods are documented for the full \pkg{sampSurf}
package in \code{\link[sampSurf:plot-methods]{plot}}.

}

\section{Methods}{

  \describe{


	\item{\strong{\code{signature(x = "ssBigBAF", y = "missing")}}}{%\cr

	  This will plot one or more sampling surfaces contained in an
	  \code{\linkS4class{ssBigBAF}} object. Optionally, it can plot a
	  set of sampled points as an overlay. }
  
  \item{}{\bold{usage\ldots}
    \preformatted{plot(x, 
     whichSS = slotNames(x)[grep('ss.', slotNames(x))],
     sampleCells = NULL,
     titles = NULL,
     namesAsTitles = TRUE,
     ... )}
   
	\itemize{
		 
	  \item{\code{x}:}{ An object that is a subclass of \dQuote{ssBigBAF}. }
		 
	  \item{\code{whichSS}:}{ The slot name(s) of the sampling surfaces
	    to plot. See \code{\linkS4class{ssBigBAF}} for the slot names. }
		 
	  \item{\code{sampleCells}:}{ A vector of cell numbers (sample
	    points) to highlight, or \code{NULL} for no points. See
	    \link{showCells} for more information. }
		 
	  \item{\code{titles}:}{ A vector (of \code{length(whichSS)}) of
	    titles that will be shown on the individual sampling surface
	    subplots; or \code{NULL} for no titles, conditional on
	    \code{namesAsTitles}.}
		 
      \item{\code{namesAsTitles}:}{ If \code{titles} is \code{NULL} or
	    not the correct length, then this kicks in as; \code{TRUE}: use
	    the slot names as titles (default); \code{FALSE}: no titles will
	    be plotted. }
	  
	  \item{\dots :}{ Other graphics arguments to be passed on to
	    the \code{\link[sampSurf:plot-methods]{plot}} method for \dQuote{sampSurf}
	    objects. }
	  
	} %itemize
 } %item plot


} %describe

} %Methods
\keyword{methods}
%\keyword{ ~~ other possible keyword(s) ~~ }
