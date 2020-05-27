\name{monte}
\alias{monte}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Generate Objects of Class "\code{\linkS4class{monteDoubleSampling}}"
}

\description{ This is the \sQuote{phantom} generic for generating
  objects of class "\code{monteDoubleSampling}." Phantom because in
  reality the generic is defined in \pkg{sampSurf} and this
  documentation just links that generic to the double sampling
  framework; i.e., there is no new \code{monte} generic definition in
  package \pkg{ssExtra}.  There are currently two methods corresponding
  to this generic for \emph{double sampling}, the main one being for
  objects of class \code{\linkS4class{ssBigBAF}}; the full documentation
  for both methods may be found in \code{\link{monte-methods}}.}

\usage{
monte(object, ...)
}

\arguments{
  
  \item{object}{ Signature argument, which differs for each method. This
    specifies the population from which samples will be drawn.}
  
  \item{\dots}{  See methods. }
}

\details{ The methods associated with this generic should be used to
  construct objects of class "\code{\linkS4class{monteDoubleSampling}}."
  These objects are specifically designed to hold information about
  Monte Carlo experiments where one resamples from a known population to
  infer efficiency and perhaps locate any bias in different sampling
  estimators. The constructor methods can be used to look at traditional
  normal theory and bootstrap confidence intervals in terms of nominal
  capture rates for the population mean.

  As already noted, the \code{monte} generic is actually defined in
  \pkg{sampSurf} (see \code{\link[sampSurf:sampSurf]{monte}} for
  details) and has been appropriated here for use in double sampling
  designs. The structure of the double sampling objects is quite
  different than those for the original \pkg{sampSurf} version. Thus,
  there is no new \code{monte} generic defined, but there are new
  methods defined for the generic in this package for double sampling
  designs. See the returned value below for the difference.
  
}

\value{
A valid object of class \dQuote{\code{\linkS4class{monteDoubleSampling}}}.
}

\references{
The \sQuote{\dQuote{monte}: When is \emph{n} Sufficiently Large?}
vignette in \pkg{sampSurf}.\cr

Please also see the User's Guide vignette for the \pkg{ssExtra} package.\cr
}

\author{
Jeffrey H. Gove %, \email{jhgove@unh.edu}
}

\note{

  Please note that there are two methods currently associated with the
  \code{monte} generic that can be used for double sampling. However,
  only the main method whose signature argument is of class
  \dQuote{ssBigBAF} should normally be called by the user. The second
  constructor method takes and object of class
  \dQuote{monteBigBAF}. This may seem a little strange, as this is the
  class that we are actually building in the first constructor. The key
  here is that the point-based Delta method (\acronym{PBDM}), or
  \dQuote{Lynch's method}, for big \acronym{BAF} sampling variance
  estimation, was not developed until after the main code was
  written. It seemed best not to add to the rather long original
  constructor method for \dQuote{ssBigBAF} objects or the code would be
  very difficult to follow, indeed much more so than it already
  is. Also, adding the code through a simple function would have
  entailed passing a large number of arguments, or a a large \code{list}
  object to the function to add the required components for the
  \acronym{PBDM}. The simplest strategy was to build an \emph{almost}
  complete object with everything but the \acronym{PBDM} variance
  components in the main (\dQuote{ssBigBAF}-based) constructor, pass
  this \emph{invalid} object to the \dQuote{monteBigBAF}-based
  constructor method, and return from that a \emph{valid} object of
  class \dQuote{monteBigBAF}.

  This process does indeed work well. Please note that calling the
  \dQuote{monteBigBAF}-based constructor with a valid object will not
  hurt anything as it will simply write over the contents of the
  \acronym{PBDM} additions. So if one is inclined to do so, one can with
  no harm done other than a little wasted \acronym{CPU} time; and there
  may be some occasion to do so, but it is largely discouraged.  The
  upshot is that one should simply use the \code{monte} constructor
  method with signature for \dQuote{ssBigBAF} objects and that is all
  that is necessary to construct a valid object of
  class\dQuote{monteBigBAF}.
  
} %note


\seealso{
\code{\link[sampSurf:sampSurf]{monte}},
\code{\linkS4class{monteDoubleSampling}},
\code{\link{monte-methods}}   
}

\examples{
#
# this example is commented-out because it consumes a lot of 
# cpu time when checking the package--just copy and paste
# it if you want to try them out...
#
# from an existing ssBigBAF object wpAll520.bb that has
# surfaces with count BAF = 5, and big BAF = 20...
#
\dontrun{
wpAll520.mc = monte(wpAll520.bb, n=c(10,25,50), start=355, mcS=50)
hist(wpAll520.mc)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
