#---------------------------------------------------------------------------
#
#   This file holds the S4 method definitions for the generation of histograms 
#   for objects of class "monteBigBAF"...
#
#   Classes and subclasses in this file include...
#
#     1. hist method for "monteBigBAF"
#
#Author...									Date: 28-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#





#================================================================================
#  1. method for class monteBigBAF...
#
setMethod('hist',
          signature(x = 'monteBigBAF'),
function(x,
         n = NA,
         xlab = '',
         col = 'gray90',
         attribute = c('vol.bb', 'ba.bb', 'ba.ct', 'vol.ct', 'tvbar'),
         ...
        )
{
#------------------------------------------------------------------------------
#
#   histogram(s) for sample size(s) n; 
#
#   The simplest thing to do is get the information for the desired variable
#   into an object of "monteNTSample" and let the hist routines in sampSurf
#   do the rest of the work.
#
#   Note that some of the slots in "monteNTSample" are either not in "monteBigBAF"
#   or they are arranged differently. Below we fudge a few slots as the
#   information that they contain does not matter for the generation of
#   histograms.
#
#   It is possible to write a more exact coercion routine that would allow
#   picking out a particular variance, for example. This could be done in the
#   future to replace the code below. Not that "setAs" is not a good fit here
#   since it has no "..." argument, and hence, we can not pick out different
#   attributes and variances, etc. The "monteBigBAF" objects contain much more
#   information that the simpler "monteNTSample" objects.
#
#   Arguments...
#     x = an argument of class "monteBigBAF"
#     n = sample size desired, NA for all available sizes
#     xlab = x-axis label
#     attribute = for the mean, (variance and st. deviation, but does not allow
#                 picking out other statistics as they are named differently
#                 and we are only concerned with histos on the mean for now.
#     ... = passed along
#
#   29-Mar-2019, JHG.
#------------------------------------------------------------------------------
#
#   cast...
#
    attribute = match.arg(attribute)
    
    ms = new('monteNTSample', 
             mcSamples = x@mcSamples,
             n = x@n,
             fpc = x@fpc,
             means = x@means[[attribute]], 
             
             #none of the variance information is used in making histograms...
             vars = x@vars[['vol.bb']],
             stDevs = x@stDevs[['vol.bb']],
             varMeans = x@varMeans[['delta']],
             stErrs = x@stErrs[['delta']],
             lowerCIs = x@lowerCIs[['delta']],
             upperCIs = x@upperCIs[['delta']],
             caught = x@caught[['delta']],
             caughtPct = x@gm.all$gm.caught['delta',],
             stats = as.data.frame(x@gm.all$gm.means),    #not correct, a kludge for now
             alpha = x@alpha,
             t.values = x@t.values,
             replace = x@replace,
             ranSeed = x@ranSeed
            )
            
      hg = hist(ms, n=n, col=col, xlab=xlab, ...)


    return(invisible(hg))

}    #hist for 'monte'
) #setMethod

