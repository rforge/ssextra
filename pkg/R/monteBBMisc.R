#---------------------------------------------------------------------------
#
#   Methods for generic show(), summary(), etc. for class monteBigBAF...
#     (1) show()
#     (2) summary()
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



#================================================================================
#  1. show() method for base class monteBigBAF...
#
setMethod('show',
          signature(object = 'monteBigBAF'),
function(object)
{
    return(summary(object, succinct=TRUE))
}   #show for 'monteBigBAF'
) #setMethod

  

#================================================================================
#  2. summary() method for class monteBigBAF...
#
setMethod('summary',
          signature(object = 'monteBigBAF'),
function(object,
         succinct = FALSE,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary...
#------------------------------------------------------------------------------
#
    if(!succinct) {
      cat('\nObject of class:', class(object))
      if(nchar(object@description) > 0) {
        .StemEnv$underLine(60)
        if(!is.na(object@description))
          cat(object@description, fill=60)
        .StemEnv$underLine(60, prologue='')
      }
      cat('Estimate attribute =', object@estimate)
      cat(' (except ba.* means)')
      if(!all(is.na(object@n))) {
        cat('\nSample sizes (n) = ')
        cat(object@n, sep=', ')
      }
      cat('\nNumber of Monte Carlo samples =', object@mcSamples)
      if(object@boot)
        cat('\nNumber of bootstrap samples =', object@numBSS, '\n')
      else
        cat('\nNo bootstrap information available.\n')
    } #if !succinct

    cat('\nMonte Carlo results...')
    cat('\n--Means...\n')
    print(object@gm.all$gm.means)
    if(!succinct) {
      cat('\n--Variances of the means...\n')
      print(object@gm.all$gm.varMeans)
      cat('\n--Standard Errors of the means...\n')
      print(object@gm.all$gm.stErrs)
    }
    cat('\n--Percent catch...\n')
    print(object@gm.all$gm.caught)

    cat('\n')
    return(invisible())
}   #summary for 'monteBigBAF'
) #setMethod
    



