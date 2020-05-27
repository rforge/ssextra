.samplingVarMeans = function(ssBB.mc)
{
#---------------------------------------------------------------------------
#
#   This little routine will calculate the approximate sampling variances
#   of the means from the MC replications, which are determined in 
#   the "monteBigBAF" constructor.
#
#   Note that this can be used as the "true" variance of the mean for
#   each sampling surface as it will converge with increasing MC reps.
#
#   In this case, to get the standard errors, just take the square root of the
#   entire sampVars df.
#
#   Arguments...
#     ssBB.mc = a "monteBigBAF" object from monte
#
#   Returns...
#     -- a data frame with the summary results
#
#**>Note: currently it excludes the tvbar, but I do not like using [1:4] here
#         as things may change, and why not include these??? Thinking...
#
#Author...									Date: Nov-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   just a quick check...
#
    if(!is(ssBB.mc, 'monteBigBAF'))
      stop('ssBB.mc object not of class "monteBigBAF"')

#
#   simpler to get these here than pass them...
#
    n = ssBB.mc@n
    n.names = names(n)
    nn = length(n)
    meanNames = names(ssBB.mc@means)          #include "tvbar" in the mix

#
#   set up the data frame for the sampling variances of the means...
#
    sampVars = data.frame(matrix(NA_real_, nrow = length(meanNames), ncol = nn))
    colnames(sampVars) = n.names
    rownames(sampVars) = meanNames
    
#
#   loop though each data frame of estimated means & calculate the vars for each n...
#
    for(s in seq_along(meanNames)) 
      sampVars[s,] = apply(ssBB.mc@means[[meanNames[s]]], 2, var, na.rm = TRUE)

#
#   wrap the variances and standard errors of the means in a list like other quantities...
#
    sm.all = list()                     
    sm.all$sm.varMeans = sampVars
    sm.all$sm.stErrs = sqrt(sampVars)
    
    return(sm.all)
}   #.samplingVarMeans


