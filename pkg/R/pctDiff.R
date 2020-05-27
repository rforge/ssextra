pctDiff = function(est,
                   true,
                   ...
                  )
{
#---------------------------------------------------------------------------
#
#   Just a simple little function that I had in-line but thought it 
#   could be more useful to have a separate copy.
#
#   Arguments...
#     est = the estimate
#     true = the "true" value
#
#   Returns...
#     simple percent difference
#
#Author...									Date: 11-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
      return((est-true)/true*100)
}   #pctDiff


