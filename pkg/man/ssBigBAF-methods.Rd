\name{ssBigBAF-methods}
\docType{methods}
\alias{ssBigBAF-methods}
\alias{ssBigBAF,sampSurf,sampSurf,sampSurf,sampSurf-method}

\title{ Methods for "ssBigBAF" Object Construction in Package \pkg{ssExtra}}

\description{ There is currently one available method for the generic
 \code{\link{ssBigBAF}}. This method generates valid objects of the associated
 "\code{\linkS4class{ssBigBAF}}" class.
}

\section{Methods}{
\describe{

\item{\strong{\code{signature(ss.bb.vol = "sampSurf", ss.bb.ba = "sampSurf", ss.ct.ba = "sampSurf", ss.ct.vol = "sampSurf")}}}{%\cr

  This is the main constructor method. There are validity checks
  performed on the individual \code{"sampSurf"} objects passed to make
  sure they are consistent in size, and store the correct
  information. These include, but are not limited to: the tract extents
  and units; all objects must share the same tree population (exactly);
  and each must be based on "\code{\linkS4class{horizontalPointIZ}}" sampling
  with matching \acronym{BAF}s for the two count surfaces, and likewise for volume.  }

  
  \item{}{\bold{usage\ldots}
    \preformatted{ssBigBAF(ss.bb.vol,
         ss.bb.ba,
         ss.ct.ba,
         ss.ct.vol,
         description = 'Big BAF sampling surface object',
         runQuiet = FALSE,
         ...)}
	   
	\itemize{
		 
	  \item{\code{ss.bb.vol}:}{ Object of class \code{"sampSurf"} that
        holds the volume surface for the large (big) \acronym{BAF} angle gauge.}

 	  \item{\code{ss.bb.ba}:}{ Object of class \code{"sampSurf"} that
        holds the basal area surface for the big \acronym{BAF} angle gauge.}

 	  \item{\code{ss.ct.ba}:}{ Object of class \code{"sampSurf"} that holds
        the basal area surface for the small/count \acronym{BAF}.}
	  
	  \item{\code{ss.ct.vol}:}{ Object of class \code{"sampSurf"} that
        holds the volume surface for the count \acronym{BAF}.}
	   
	  \item{\code{description}:}{ A character description for the return
	  object.}
		 
	  \item{\code{runQuiet}:}{ \code{TRUE}: no feedback; \code{FALSE}:
	    feedback if available. }
		
	  \item{... :}{ Other arguments to be passed along perhaps
	  eventually, right now they are \dQuote{gobbled.} }

	} %itemize
 } %item

} %describe
} %section
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
