\name{ssBigBAF-class}
\Rdversion{1.1}
\docType{class}
\alias{ssBigBAF-class}
%% \alias{plot,ssBigBAF,missing-method}
\alias{show,ssBigBAF-method}
\alias{summary,ssBigBAF-method}

\title{Class \code{"ssBigBAF"}}

\description{
  
  This class defines a subclass of \code{"\linkS4class{ssDoubleSampling}"} for use with
  so-called \dQuote{Big \acronym{BAF}} sampling.

}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("ssBigBAF", ...)}.  However, it is
  simpler and less error-prone to use the constructor of the same name:
  \code{\link{ssBigBAF}}.

}



\section{Slots}{
  
  In addition to those slots that are defined within the
  \code{"\linkS4class{ssDoubleSampling}"}, the following are
  available\ldots
  
  \describe{
	
    \item{\code{ss.bb.vol}:}{Object of class \code{"sampSurf"} that
    holds the volume surface for the large (big) \acronym{BAF} angle gauge. }
	
    \item{\code{ss.bb.ba}:}{Object of class \code{"sampSurf"} that holds
    the basal area surface for the big \acronym{BAF} angle gauge.}
	
    \item{\code{ss.ct.ba}:}{Object of class \code{"sampSurf"} that holds
    the basal area surface for the small/count \acronym{BAF}.}
	
    \item{\code{ss.ct.vol}:}{Object of class \code{"sampSurf"} that
      holds the volume surface for the count \acronym{BAF}. }

    \item{\code{csl.bb}:}{Object of class \code{"ssCellStemList"} 
       corresponding to the big \acronym{BAF} surfaces.	}

    \item{\code{csl.ct}:}{Object of class \code{"ssCellStemList"} 
      corresponding to the count \acronym{BAF} surfaces.	}

  } %describe
} %section

\section{Extends}{
Class \code{"\linkS4class{ssDoubleSampling}"}, directly.
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "ssBigBAF", y = "missing")}: plot
      the surface(s). }
	
    \item{show}{\code{signature(object = "ssBigBAF")}: used to print the
    object. }
	
    \item{summary}{\code{signature(object = "ssBigBAF")}: for an object summary. }
	 }
}
   
%% \references{
%% %%  ~~put references to the literature/web site here~~
%% }

\author{
  
\packageAuthor{ssExtra}

Maintainer: \packageMaintainer{ssExtra}

}

%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  The parent (super) class: \code{\linkS4class{ssDoubleSampling}}
}

\examples{
showClass("ssBigBAF")
}
\keyword{classes}
