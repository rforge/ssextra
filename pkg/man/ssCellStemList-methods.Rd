\name{ssCellStemList-methods}
\docType{methods}
\alias{ssCellStemList-methods}
\alias{ssCellStemList,sampSurf-method}

\title{ Methods for "ssCellStemList" Object Construction in Package
  \pkg{ssExtra} }

\description{

There is currently one available method for the generic
 \code{\link{ssCellStemList}}. This method generates valid objects of
 the associated "\code{\linkS4class{ssCellStemList}}" class.

}

\section{Methods}{
  
\describe{

\item{\strong{\code{signature(ss = "sampSurf")}}}{


  This is the main constructor method. }

  
  \item{}{\bold{usage\ldots}
    \preformatted{ssCellStemList(ss,
	       description = 'list of cells in stems & stems in cells',
               ...)}

	   
	\itemize{
		 
	  \item{\code{ss}:}{ An object of class \code{"sampSurf"}. }
	   
	  \item{\code{description}:}{ A character description for the return
	  object.}
		
	  \item{... :}{ Other arguments to be passed along perhaps
	  eventually, right now they are \dQuote{gobbled.} }

	} %itemize
 } %item


} %describe
} %section

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
