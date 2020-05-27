\name{ssBigBAF}
\alias{ssBigBAF}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Generate Objects of Class "\code{\linkS4class{ssBigBAF}}"
}

\description{ This generic function has one method based on the four
signature formal arguments. It is used as a constructor function for
objects that are of class "\code{\linkS4class{ssBigBAF}}" and should be
used in preference to \code{\link{new}} to insure a valid object.  }

\usage{
ssBigBAF(ss.bb.vol,
         ss.bb.ba,
         ss.ct.ba,
         ss.ct.vol,
         ...)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{ss.bb.vol}{ Object of class \code{"sampSurf"} that
    holds the volume surface for the large (big) \acronym{BAF} angle gauge.
  }
  
  \item{ss.bb.ba}{Object of class \code{"sampSurf"} that holds
    the basal area surface for the big \acronym{BAF} angle gauge.
  }

  \item{ss.ct.ba}{Object of class \code{"sampSurf"} that holds
    the basal area surface for the small/count \acronym{BAF}.
  }
  
  \item{ss.ct.vol}{Object of class \code{"sampSurf"} that
    holds the volume surface for the count \acronym{BAF}.
  }
  
  \item{\dots}{ Formal arguments that are different for each method, see
	those for details.
  }
  
} %arguments

\details{ \code{ssBigBAF} is defined completely with respect to the
signature of its method \code{\link{ssBigBAF-methods}}.
}

\value{
An object of class "\code{\linkS4class{ssBigBAF}}".
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Jeffrey H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[sampSurf:sampSurf]{sampSurf}}
}

%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% ## The function is currently defined as
%% "ssBigBAF"
%% }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
