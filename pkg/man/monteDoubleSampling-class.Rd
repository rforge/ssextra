\name{monteDoubleSampling-class}
\Rdversion{1.1}
\docType{class}
\alias{monteDoubleSampling-class}

\title{Class \code{"monteDoubleSampling"}}

\description{
A virtual base class that allows creation of subclasses related to Monte
Carlo sampling from subclasses of the
\code{"\linkS4class{ssDoubleSampling}"} virtual class.
}

\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}

\section{Slots}{
  
  \describe{
	
    \item{\code{description}:}{Object of class \code{"character"}: A
       description of the object if desired (defaults are given for each
       subclass). }
	
    \item{\code{estimate}:}{Object of class \code{"character"}: The
      type of estimate or attribute to which Monte Carlo subsampling applies. }
  }
  
}

\section{Methods}{
No methods defined with class \dQuote{monteDoubleSampling} in the signature.
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

\seealso{
  The following subclasses and related: \code{\linkS4class{monteBigBAF}},
}

\examples{
showClass("monteDoubleSampling")
}
\keyword{classes}
