.grandMeans = function(x, n.names)
{
#---------------------------------------------------------------------------
#
#   This routine will calculate the grand means over each sample for each n.
#   The results are akin to the normal theory summary statistics; e.g., the
#   sampling variance of the means calcluated by .samplingVarMeans, but
#   are not, of course the same. They should probably converge as the number
#   of MC draws and sample size gets large.
#
#   Note that it also can contain jackknife & bootstrap results, depending
#   on what was chosen in the "monteBigBAF" object constructor.
#
#   This little function does the summary stat work--note that sapply will
#   happily drop from matrix to numeric if there is only one sample size, 
#   so we catch and correct this below.
#
#   Arguments...
#     x = a list from the "monteBigBAF" constructor for one set of estimate
#         data frames
#     n.names = the names for sample sizes
#
#   Returns...
#     -- a matrix with the summary "grand" means
#
#Author...									Date: 25-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   sapply can truncate the dimensions as noted above...
#
    res = sapply(x, colMeans, na.rm = TRUE)
    
#
#   check for dim reduction if only one sample size & restore to matrix...
#
    if(!is(res, 'matrix')) {        #if TRUE, it's a numeric vector
      res = as.matrix(res)
      rownames(res) = names(x)
      colnames(res) = n.names
    }
    else
      res = t(res)                  #to fix the way sapply returns the results
      
    return(res)
}   #.grandMeans

