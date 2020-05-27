#---------------------------------------------------------------------------
#
#   Methods for generic show(), summary(), etc. for class ssBigBAF...
#     (1) show()
#     (2) summary()
#
#Author...									Date: 22-Mar-2019
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
#  (1) show() method for class "ssBigBAF"...
#
setMethod('show',
          signature(object = 'ssBigBAF'),
function(object)
{
    return(summary(object))

    cat('\n')
    return(invisible())
}   #show for 'ssBigBAF'
) #setMethod




#================================================================================
#  (2) summary method for class ssBigBAF...
#
setMethod('summary',
          signature(object = 'ssBigBAF'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple summary of common items...
#------------------------------------------------------------------------------
#
    cat('\nObject of class:', class(object))
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    #some simple little extraction functions useful in the report function...
    theBAF = function(ss) ss@izContainer@iZones[[1]]@angleGauge@baf
    theUnits = function(ss) ss@izContainer@iZones[[1]]@angleGauge@units
    theIZs = function(ss) class(ss@izContainer@iZones[[1]])

    ssNames = slotNames(object)[grep('ss.', slotNames(object))]
    #prints a line of the report for a given sampling surface object...
    report = function(ss, name, bc) {
               cat('\n  ', name, ': ', bc, ' BAF ', 
                           ss@estimate, ' surface, ', 
                           theIZs(ss), ', ', 
                           theBAF(ss), ' BAF (', theUnits(ss), ')', 
                           sep=''
                  )
               return(invisible())
    } #report
    
    cat('\nSampling surface slot contents/estimates...')
    nSS = length(ssNames)
    for(i in seq_len(nSS)) {
      bb.ct = ifelse(length(grep('bb', ssNames[i])) > 0, 'Big', 'Count')
      ss = parse(text = paste0('object@', ssNames[i]))
      report(eval(ss), ssNames[i], bb.ct)
    }

    cat('\nNumber of stems =', length(object@ss.bb.vol@izContainer@iZones) )

#
#   the following is based on the mean surface basal areas, and takes into
#   consideration multiple trees per point because when divided by the appropriate
#   baf, this gives the total tree count for the entire surface, and thus the 
#   expected ratio; it will be ever so slightly different than baf.v/baf.c...
#
    baf.c = object@ss.ct.vol@izContainer@iZones[[1]]@angleGauge@baf
    baf.v = object@ss.bb.vol@izContainer@iZones[[1]]@angleGauge@baf
    n.bb = object@ss.bb.ba@surfStats$mean/baf.v
    n.ct = object@ss.ct.ba@surfStats$mean/baf.c
    ss.ratio = n.ct/n.bb
    cat('\nMean BA surface-based sampling ratio count:bb =', ss.ratio)
    baf.ratio = baf.v/baf.c
    cat('\nTrue baf-based sampling ratio count:bb =', baf.ratio)

    cat('\n\nPlease use "summary" on individual components (slots) for more details.')

    cat('\n')
    
    return(invisible())
}   #summary for 'ssBigBAF'
) #setMethod


