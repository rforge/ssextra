\name{ssExtra-package}
\alias{ssExtra-package}
%\alias{ssExtra}
\docType{package}

%
% note that the \package*{} commands below extract the relevant
% information from the package DESCRIPTION file. So build the rest
% around these or delete them...

\title{
  
\packageTitle{ssExtra}: Double Sampling, &c.

}

\description{   %\packageDescription{ssExtra}

  This package provides some potentially useful helper functions for
  working with objects created from the \pkg{sampSurf} package. It also
  provides some S4 classes that extend the use of sampSurf to, e.g.,
  double sampling.  Currently, the main double sampling application is
  the Monte Carlo addition to the generic \code{\link[sampSurf]{monte}}
  constructor for \dQuote{Big \acronym{BAF}} sampling.

  At this point the \pkg{sampSurf} package is quite large. Rather than
  continually adding to it until is capsizes under its own weight,
  perhaps the better option is to extend it with other smaller, more
  managable packages. This was done with the wavelet additions in the
  \pkg{ssWavelets} package, and thus is also used here. The routines in
  both \pkg{ssWavelets} and \pkg{ssExtra} are more specialized in
  general (though there are some interesting helper routines too) and
  thus placing them in other packages to be loaded only as needed seems
  to be a reasonable approach.

}

\details{
  
\tabular{ll}{
Package: \tab ssExtra\cr
Type: \tab Package\cr
Version: \tab 0.1-2\cr
Date: \tab 12-May-2020\cr
License: \tab GPL (>= 3)\cr
LazyLoad: \tab yes\cr
}

Please see the package vignette for some examples of how to use the
routines in this package in conjunction with those in \pkg{sampSurf}.

} %details

\section{Classes for Use with \pkg{sampSurf}\ldots}{

  The following S4 classes are defined within the \pkg{ssExtra} package
  to augment those in \pkg{sampSurf} with extended functionality.
  Objects of the non-virtual classes below can be created using
  constructor functions of the same name, see the \emph{Object
  Constructors} section for details.

%
  
  \subsection{The \dQuote{ssBigBAF} class\ldots}{ This is a simple class
    made up of a related set of \dQuote{\code{\linkS4class{sampSurf}}}
    objects.

	\tabular{ll}{
      \code{\linkS4class{ssDoubleSampling}} \tab Virtual base class for the following\cr
      \code{\linkS4class{ssBigBAF}} \tab Class structure for \dQuote{Big
      \acronym{BAF}} sampling\cr
    }
  
  }%ssBigBAF

%
  
  \subsection{The \dQuote{ssCellStemList} class\ldots}{ This is a simple class
    with no subclasses.

	\tabular{ll}{
      \code{\linkS4class{ssCellStemList}} \tab Lists of Cells-in-Stems
      (inclusion zones) or Stems-in-Cells.
    }

	Objects of this class are used in slots within
	\code{"\linkS4class{ssBigBAF}"} objects; one slot for each of the
	big \acronym{BAF} and count surfaces. It is unlikely, but not
	impossible that this class might be of use in casual application without the
	link to big \acronym{BAF} sampling.
	
  }%ssCellStemList

  
  \subsection{The \dQuote{monte} classes\ldots}{ These classes are more
    complex extensions of the \pkg{sampSurf} \code{\linkS4class{monte}}
    class. These new classes are used for Monte Carlo subsampling from
    objects of a subclass of \dQuote{ssDoubleSampling}.

	\tabular{ll}{
      \code{\linkS4class{monteDoubleSampling}} \tab Virtual base class for the following\cr
      \code{\linkS4class{monteBigBAF}} \tab Class structure for \dQuote{Big
      \acronym{BAF}} sampling\cr
    }

   Note that these classes are \emph{not} directly related though
   inheritance to those in the \dQuote{monte} class structure within the
   standard \pkg{sampSurf} package. Because we are interested in several
   different quantities such as basal area and volume on large and small
   \acronym{BAF} sampling surfaces, and the comparison of different
   variance estimators, much of the structure of the current
   \dQuote{monte}-type objects is in the form of lists of data frames,
   rather than simple data frames as in the original \pkg{sampSurf}
   implementation. Hence the use of inheritance was not
   possible. However, it is possible to manually coerce the more complex
   objects to the simpler structure (with some necessary loss of
   information); see, for example, the \code{hist} method below.

   In additon, it is possible that some of the slots in the
   "\code{monteBigBAF}" subclass might migrate to the virtual
   "\code{monteDoubleSampling}"superclass in the future. This will not
   change anything for the user. The plan for future
   "\code{monteDoubleSampling}" subclasses was not really considered at
   the time of the coding for the current classes, therefore much of
   what might end up being common slot-based functionality could be
   moved eventually.
	
  }%monteDoubleSampling

} %S4 classes









%--------------------------------------------------------------------------
%  constructors section...
%
\section{Object Constructors}{
  
  For each of the classes defined in the section above, we must be able
  to create objects that can be used in \R. This is done using
  class-specific \dQuote{constructor} (sometimes called
  \dQuote{generator}) methods that take the drudgery away from creating
  what can be somewhat complicated new object instances. There may be
  more than one constructor for a given class of object, and these are
  differentiated by the method signature; see the links provided below
  for more details\ldots

%
  \subsection{\dQuote{ssBigBAF} class constructor\ldots}{

  \tabular{ll}{
	\code{\link{ssBigBAF}} \tab Constructor for individual
	  \dQuote{ssBigBAF} objects\cr   
  }

  }%ssBigBAF

%
  \subsection{\dQuote{ssCellStemList} class constructor\ldots}{

  \tabular{ll}{
	\code{\link{ssCellStemList}} \tab Constructor for individual
	  \dQuote{ssCellStemList} objects\cr   
  }

  }%ssCellStemList

%
  \subsection{\dQuote{monteBigBAF} class constructor and helper functions\ldots}{

  \tabular{ll}{
	\code{\link{monte}} \tab Constructor for individual
	\dQuote{monteDoubleSampling} subclass objects\cr
	\code{\link{monteStatsBB}} \tab Calculate expanded point-wise statistics for
	individual Monte Carlo sample draws\cr
	\code{\link{monteTreeStatsBB}} \tab Calculate unexpanded tree-wise statistics for
	individual Monte Carlo sample draws\cr
	\code{\link{monteBootBB}} \tab Calculate jackknife and bootstrap
	samples on individual Monte Carlo replicates\cr
  }

  }%monteBigBAF

	
} %constructors






%--------------------------------------------------------------------------
%  summary & plotting section...
%

\section{Summary and plotting methods}{
  
  The objects created from the above classes have graphical
  content made possible by using classes from the \code{\link{sp}} or
  \code{\link{raster}} packages. Therefore, methods have been added to
  the \code{\link{plot}} generics to allow for graphical display of
  objects. In addition, \code{\link{summary}} and \code{\link{show}}
  methods are also available (sometimes producing the same result) for
  printing information within the object\ldots

  \tabular{ll}{
	\code{\link[ssExtra:hist-methods]{hist}} \tab Create a histogram of \dQuote{monteBigBAF} results\cr	
	\code{\link[ssExtra:plot-methods]{plot}} \tab  To plot, e.g., a \dQuote{ssBigBAF} object\cr
	\code{\link{show}} \tab  To succinctly print one of the objects\cr
	\code{\link{summary}} \tab  To summarize one of the objects\cr
	\code{\link{showCells}} \tab Can be used with a number of different
	objects to plot a sample of grid cells (sample points)\cr
  } %tabular
  
} %summary and plotting



%--------------------------------------------------------------------------
%  Other methods section...
%


\section{Miscellaneous methods}{

  The following list a few methods that may be of some help in working
  with \pkg{sampSurf} and \pkg{ssExtra}\ldots

  \tabular{ll}{
	\code{\link{initTract}} \tab Creates a simple
	    \code{\linkS4class{bufferedTract}} object for simulations\cr
	\code{\link{drawTreePop}} \tab Creates a reasonable synthetic
	    northern hardwoods tree population\cr
	\code{\link{makePop}} \tab Creates a basic population for
        simulations on a \dQuote{bufferedTract} using the above\cr
    \code{\link{getTreesBB}} \tab Coerces a
	  \code{\linkS4class{standingTrees}} object to a data frame while adding
	  tree VBARs and basal area \cr
    \code{\link{listCellStem}} \tab Handy routine to list the the
	  sampling surface cell numbers within each tree's inclusion zone, and
	  vice versa \cr
	\code{\link{testBoot}} \tab Compare different bootstrap confidence
      intervals \cr
	\code{\link{compareBoot}} \tab Compare sets of bootstrap confidence
      intervals \cr
	\code{\link{createBBNH}} \tab Create an \dQuote{ssBigBAF} object for
	the synthetic northern hardwood population \cr
	\code{\link{pctDiff}} \tab Simple percent difference function \cr
	\code{\link{ratioVariance}} \tab Used to calculate the ratio
	variance for the VBAR ratio \cr
  } %tabular
  
} %miscellaneous methods


%----------------------------------------------------------------------------------
\author{
  
\packageAuthor{ssExtra}

Maintainer: \packageMaintainer{ssExtra}
}

\references{

\describe{
%

\item{}{Gove, J. H., Gregoire,T. G., Ducey, M. J., and Lynch,
   T. B. 2020. A Note on the Estimation of Variance for Big BAF
   Sampling. \emph{Forest Ecosystems}, Submitted.}

}%describe
  
}%references

\keyword{ package }

\seealso{
  
 Please see the \code{\link[sampSurf:sampSurf-package]{sampSurf}}
 package as well as the \code{\link{sp}} and
 \code{\link[raster:raster-package]{raster}} packages for details on the
 underlying grid and polygon classes and methods.

}

%% \examples{
%% ~~ simple examples of the most important functions ~~
%%}
