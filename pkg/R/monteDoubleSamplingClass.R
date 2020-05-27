#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the "monteDoubleSampling" 
#   virtual and subclass(es)
#
#   Classes and subclasses in this file include...
#
#     1. monteDoubleSampling means --  VIRTUAL
#     2. monteBigBAF          --  subclass of monteDoubleSampling
#
#   I probably should have set up at least a pair of SS slots in the virtual
#   class. The problem is, what to call them. It seems like the application
#   will vary based on what X & Y are attribute-wise in any one situation
#   so should they just be called ss.X and xx.Y, or have some subclass
#   fill in the details with a more appropriate set of names? Since this
#   was a means to an end w/r to Big BAF sampling, I decided to just leave
#   the virtual base class essentially empty, and define slots with meaningful
#   names in the subclasses for now. This could be changed at the later
#   dat, but probably will not be.
#
#Author...									Date: 27-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#



#=================================================================================================
#
#  1. define the monteDoubleSampling class...
#
setClass('monteDoubleSampling',
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',
                   estimate = 'character'                   #one of .StemEnv$puaEstimates
                   ),
                   
    contains = 'VIRTUAL',                                      #note!
    
    #some defaults for validity checking...
    prototype = list(description = 'Monte Carlo Double Sampling',
                     estimate = character(0)
                    ),
                    
    validity = function(object) {

#              no checks for now...
                   
                 return(TRUE)
               } #validity check
) #class monteDoubleSampling 




##setClassUnion('listOrNULL', c('list', 'NULL')) #see below

#=================================================================================================
#
#  2. the monteBigBAF class is a direct descendant/subclass of 'monteDoubleSampling'...
#
#     What follows is very terse w/r to explanation of what the different "list" slots contain.
#     They are not all standardized with the same content. Please see the constructor for more
#     information on what each contains.
#
#     The first several slots are the same as what appears in the "monteSample" virtual class
#     definition. These could (and should) probably be moved up into the virtual superclass
#     "monteDoubleSampling".
#
#     The other slots are similar to "monteSample" as well, but due to the fact that we have four
#     sampling surfaces and several variance estimators, the information is in lists. However,
#     the data frames within the lists all have standard form similar to those in "monteSample"
#     and its subclasses.
#
#     Note that if 'boot=TRUE' then the lists will contain estimates from the jackknife and
#     bootstrap as well; otherwise, this information will be missing.
#

setClass('monteBigBAF',
#
#  slots for the class and its subclasses...
#
    representation(mcSamples = 'numeric',               #number of Monte Carlo samples
                   n = 'numeric',                       #vector of sample sizes
                   fpc = 'numeric',                     #finite pop correction by sample size
                   alpha = 'numeric',                   #two-tailed alpha level
                   replace = 'logical',                 #TRUE: MC sample with replacement
                   ranSeed = 'numeric',                 #initial random seed value
                   t.values = 'numeric',                #two-tailed t-value for each level m
                   boot = 'logical',                    #TRUE: jack & boot estimates included
                   numBSS = 'numeric',                  #number of boostrap samples
                   
                   #the following lists contain data frames with dim: (mcSamples x length(n))
                   means = 'list',                      #sample means by sample size
                   vars = 'list',                       #sample variances by sample size
                   stDevs = 'list',                     #sample standard deviations by sample size
                   varMeans = 'list',                   #variances of mean estimates
                   stErrs = 'list',                     #standard errors of the means
                   lowerCIs = 'list',                   #lower limit of the confidence intervals
                   upperCIs = 'list',                   #upper limit of the confidence intervals
                   caught = 'list',                     #logical data frames of caught stats
                   otherVarms = 'list',                 #variance of mean for other quantities
                   n.tvbar = 'list',                    #number of vbar trees, for count & BB
                   corrs = 'list',                      #correlations for BB BA and tree vbars
                   covs = 'list',                       #covariances for BB BA and tree vbars
                   
                   gm.all = 'list',                     #grand-summary (means) of above
                   sm.all = 'list',                     #MC sampling variance of the means
                   mc.samples = 'list'                  #n-list of (n x mcSamples) data frame
                   
                  ),
                   
                   
    contains = 'monteDoubleSampling',                   #descendent of monteDoubleSampling
    
    #some defaults for validity checking...
    prototype = list(description = 'Big BAF Sampling Monte Carlo Results',
                     estimate = .StemEnv$puaEstimates$volume,
                     mcSamples = NA_real_,
                     n = NA_real_,
                     alpha = NA_real_,
                     fpc = NA_real_,
                     replace = NA,
                     ranSeed = NA_real_,
                     t.values = NA_real_,
                     boot = NA,
                     numBSS = NA_real_,
                     #
                     means = list(),
                     vars = list(),
                     stDevs = list(),
                     varMeans = list(),
                     stErrs = list(),
                     lowerCIs = list(),
                     upperCIs = list(),
                     caught = list(),
                     otherVarms = list(),
                     n.tvbar = list(),
                     corrs = list(),
                     covs = list(),
                     gm.all = list(),
                     sm.all = list(),
                     mc.samples = list()
                    ),
 
                    
    validity = function(object) {
    
          #check here that all lists contain data frames (gm.all are matrices)...
          if(!all(sapply(object@means, is, 'data.frame')))
            stop('***>All components in the \'means\' list must be of class "data.frame"')
          if(!all(sapply(object@vars, is, 'data.frame')))
            stop('***>All components in the \'vars\' list must be of class "data.frame"')
          if(!all(sapply(object@stDevs, is, 'data.frame')))
            stop('***>All components in the \'stDevs\' list must be of class "data.frame"')
          if(!all(sapply(object@varMeans, is, 'data.frame')))
            stop('***>All components in the \'varMeans\' list must be of class "data.frame"')
          if(!all(sapply(object@stErrs, is, 'data.frame')))
            stop('***>All components in the \'stDevs\' list must be of class "data.frame"')
          if(!all(sapply(object@lowerCIs, is, 'data.frame')))
            stop('***>All components in the \'lowerCIs\' list must be of class "data.frame"')
          if(!all(sapply(object@upperCIs, is, 'data.frame')))
            stop('***>All components in the \'upperCIs\' list must be of class "data.frame"')
          if(!all(sapply(object@caught, is, 'data.frame')))
            stop('***>All components in the \'caught\' list must be of class "data.frame"')
          if(!all(sapply(object@otherVarms, is, 'data.frame')))
            stop('***>All components in the \'otherVarms\' list must be of class "data.frame"')
          if(!all(sapply(object@n.tvbar, is, 'data.frame')))
            stop('***>All components in the \'n.tvbar\' list must be of class "data.frame"')
          if(!all(sapply(object@corrs, is, 'data.frame')))
            stop('***>All components in the \'corrs\' list must be of class "data.frame"')
          if(!all(sapply(object@covs, is, 'data.frame')))
            stop('***>All components in the \'covs\' list must be of class "data.frame"')

          if(!all(sapply(object@gm.all, is, 'matrix')))
            stop('***>All components in the \'gm.all\' list must be of class "matrix"')  #matrix
          if(!all(sapply(object@sm.all, is, 'data.frame')))
            stop('***>All components in the \'sm.all\' list must be of class "data.frame"')
          if(!all(sapply(object@mc.samples, is, 'data.frame')))
            stop('***>All components in the \'mc.samples\' list must be of class "data.frame"')
                   
      return(TRUE)
    } #validity check
    
) #class monteBigBAF


