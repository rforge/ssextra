#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   ssBigBAF class...
#
#   The methods include signatures for...
#     Signature: "object", "object", "object", "object"
#     1. "sampSurf", "sampSurf", "sampSurf", "sampSurf": Takes four sampSurf
#         objects
#
#
#Author...									Date: 21-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#   generic definition...
#
setGeneric('ssBigBAF',  
           function(ss.bb.vol, ss.bb.ba, ss.ct.ba, ss.ct.vol, ...) standardGeneric('ssBigBAF'),
             signature = c('ss.bb.vol', 'ss.bb.ba', 'ss.ct.ba', 'ss.ct.vol')
            )




          
#================================================================================
#
# 1. Takes a set of "sampSurf" objects...
#
setMethod('ssBigBAF',
          signature(ss.bb.vol = 'sampSurf', ss.bb.ba = 'sampSurf', 
                    ss.ct.ba = 'sampSurf',  ss.ct.vol = 'sampSurf'), 
function(ss.bb.vol, 
         ss.bb.ba,
         ss.ct.ba,
         ss.ct.vol,
         description = 'Big BAF sampling surface object',
         runQuiet = FALSE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   This is the main way to create a simple "ssBigBAF" object...
#
#---------------------------------------------------------------------------
#
#   create the ssCellStemList objects...
#
    if(!runQuiet)
      cat('\nCreating the big BAF surface cell-stem IZ list...')
    csl.bb = ssCellStemList(ss.bb.vol)
    
    if(!runQuiet)
      cat('\nCreating the count surface cell-stem IZ list...')
    csl.ct = ssCellStemList(ss.ct.ba)

#
#   create the object...
#
    ssBB = new('ssBigBAF',
               description = description,
               ss.bb.vol = ss.bb.vol,
               ss.bb.ba = ss.bb.ba,
               ss.ct.ba = ss.ct.ba,
               ss.ct.vol = ss.ct.vol,
               csl.bb = csl.bb,
               csl.ct = csl.ct
              )

    if(!runQuiet)
      cat('\n')
      
    return(ssBB)
}   #ssBigBAF for 4*"sampSurf"
)   #setMethod
        


