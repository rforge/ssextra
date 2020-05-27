\name{ssDoubleSampling-class}
\Rdversion{1.1}
\docType{class}
\alias{ssDoubleSampling-class}

\title{Class \code{"ssDoubleSampling"}}

\description{
  This is a virtual class that establishes a very bare framework for
  creating subclasses that are based on double sampling. Otherewise, it
  has no utility and only one slot is inherited.
}

\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}

\section{Slots}{
  \describe{
    \item{\code{description}:}{Object of class \code{"character"}: A
       description of the object if desired (defaults are given for each
       subclass). }
  }
}

\section{Methods}{
No methods defined with class "ssDoubleSampling" in the signature.
}

%% \references{
%% %%  ~~put references to the literature/web site here~~
%% }

\author{
  
\packageAuthor{ssExtra}

Maintainer: \packageMaintainer{ssExtra}

}

\seealso{
  The following subclasses and related: \code{\linkS4class{ssBigBAF}},
}

\examples{
showClass("ssDoubleSampling")
}
\keyword{classes}
