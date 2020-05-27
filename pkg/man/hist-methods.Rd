\name{hist-methods}
\docType{methods}
%%\alias{hist,monte-method}
\alias{hist,monteBigBAF-method}

\title{ Methods for Function \code{hist} in Package \pkg{ssExtra} }

\description{

  Simple methods to generate histograms for select objects in \pkg{ssExtra}.
  
}

\section{Methods}{

  \describe{


   \item{\strong{\code{signature(x = "monteBigBAF")}}}{%\cr

	 This will generate as set of histograms for objects of class
	 \code{\linkS4class{monteBigBAF}}. One can choose the sample sizes,
	 or display subplots for all sample sizes by default.   }

     
  \item{}{\bold{usage\ldots}
    \preformatted{hist(x, 
     n = NA,
     xlab = '',
     col = 'gray90',
     attribute = c('vol.bb', 'ba.bb', 'ba.ct', 'vol.ct', 'tvbar'),
     ... )}
   
	\itemize{
		 
	  \item{\code{x}:}{ An object that is a subclass of \dQuote{monteBigBAF}. }
		 
	  \item{\code{n}:}{ A vector of sample sizes specifying the
	    histogram subfigures; \code{NA} plots all available sample
	    sizes. }
		 
	  \item{\code{xlab}:}{ The label for the x-axis. }
		 
	  \item{\code{col}:}{ A color for the histogram bars. }

	  \item{\code{attribute}:}{ The attribute means to plot: \sQuote{vol.bb} is volume
          from the Big \acronym{BAF} surface; \sQuote{ba.bb} is basal area from that
          surface; \sQuote{ba.ct} is basal area from the count surface; \sQuote{vol.ct} is
          volume from the count surface; and \sQuote{tvbar} are the raw tree VBARs. }
		
	  \item{\dots :}{ Other arguments to be passed on to the \code{hist}
	    function for objects of class \code{\linkS4class{monteNTSample}}
	    in \dQuote{sampSurf}. }
	  
	} %itemize
 } %item hist
   


}%describe
}%Methods

%% \keyword{methods}
%% \keyword{ ~~ other possible keyword(s) ~~ }
