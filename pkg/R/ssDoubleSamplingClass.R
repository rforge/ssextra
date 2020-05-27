#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the "ssDoubleSampling" 
#   virtual and subclass(es)
#
#   Classes and subclasses in this file include...
#
#     1. ssDoubleSampling  --  VIRTUAL
#     2. ssBigBAF          --  subclass of ssDoubleSampling
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
#




#=================================================================================================
#
#  1. define the ssDoubleSampling class...
#
#     Note that I could have defined four generic sampling surface slots here, and may in the
#     future. But since each application could involve different attributes, it becomes
#     potentially confusing to find good names. Perhaps something like...
#       ss.1a, ss.1b, ss.2a, ss.2b, where "1" is the primary and "2" the sub-sample; and
#     "a" and "b" stand for attributes. Seems cumbersome.
#
#     One could also make a list of four sampSurf objects, then name that list in each
#     child application corresponding to what is in them. However, this then loses the
#     type-checking structure in S4 of having them directly checked on invocation of
#     "new".
#
setClass('ssDoubleSampling',
#
#  slots for the class and its subclasses...
#
    representation(description = 'character'
                   ),
                   
    contains = 'VIRTUAL',                                      #note!
    
    #some defaults for validity checking...
    prototype = list(description = 'sampSurf Double Sampling'
                    ),
                    
    validity = function(object) {

#              no checks for now...
                   
                 return(TRUE)
               } #validity check
) #class ssDoubleSampling 





#=================================================================================================
#
#  2. the ssBigBAF class is a direct descendant/subclass of 'ssDoubleSampling'...
#
#     Note the validity check below on the classes of the SS objects. The check for 
#     "horizontalPointIZ" allows for its subclasses to work seemlessly as well. For
#     example, "criticalHeightIZ" would work for all four surfaces. Also, one could
#     conceivably mix them so that the Count BAF surfaces were "horizontalPointIZ"
#     while the Big BAF surfaces were "criticalHeightIZ" or some antithetic subclass.
#     The latter is allowed since, in practice, we would not have volumes on the
#     count plots; but these are used in the monte() simulations for checking H-T
#     variance, for example.
#
#
setClass('ssBigBAF',
#
#  slots for the class and its subclasses...
#
    representation(ss.bb.vol = 'sampSurf',                   #big BAF volume sampling surface
                   ss.bb.ba = 'sampSurf',                    #big BAF basal area sampling surface
                   ss.ct.ba = 'sampSurf',                    #count BAF basal area sampling surface
                   ss.ct.vol = 'sampSurf',                   #count BAF volume sampling surface
                   #only need one of each...
                   csl.bb = 'ssCellStemList',                #for the two big BAF surfaces
                   csl.ct = 'ssCellStemList'                 #for the two count surfaces
                  ),
                   
    contains = 'ssDoubleSampling',                           #descendent of ssDoubleSampling
    
    #some defaults for validity checking...
    prototype = list(description = 'Big BAF Sampling',
                     ss.bb.vol =  new('sampSurf'),           #dummy, generates an invalid object
                     ss.bb.ba =  new('sampSurf'),            #dummy, generates an invalid object
                     ss.ct.ba =  new('sampSurf'),            #dummy, generates an invalid object
                     ss.ct.vol =  new('sampSurf'),            #dummy, generates an invalid object
                     #more dummy objects...
                     csl.bb = new('ssCellStemList'),
                     csl.ct = new('ssCellStemList')
                    ),

                    
    validity = function(object) {
        #first, the units must be the same in all four tracts...
        if(object@ss.bb.vol@tract@units != object@ss.bb.ba@tract@units ||
           object@ss.bb.vol@tract@units != object@ss.ct.ba@tract@units ||
           object@ss.bb.vol@tract@units != object@ss.ct.vol@tract@units)
          stop('***>All tracts must have the same units of measure!')
        
        #make sure the tract extents are the same for all...
        if(!identical(extent(object@ss.bb.vol@tract), extent(object@ss.bb.ba@tract)) || 
           !identical(extent(object@ss.bb.vol@tract), extent(object@ss.ct.ba@tract)) ||
           !identical(extent(object@ss.bb.vol@tract), extent(object@ss.ct.vol@tract)))
         stop('***>Tract extents must be the same for each SS object!')

        #make sure all objects are either 'horizontalPointIZ' or a subclass, which
        #would include chs, for example...
        if(!is(object@ss.bb.vol@izContainer@iZones[[1]], 'horizontalPointIZ') ||
           !is(object@ss.bb.ba@izContainer@iZones[[1]], 'horizontalPointIZ') ||
           !is(object@ss.ct.ba@izContainer@iZones[[1]], 'horizontalPointIZ') ||
           !is(object@ss.ct.vol@izContainer@iZones[[1]], 'horizontalPointIZ') )
          stop('***>At least one sampling surface is not of class "horizontalPointIZ"!')

         
        #check for identical angleGauge objects: very strict...
        if(!identical(object@ss.bb.vol@izContainer@iZones[[1]]@angleGauge,
                      object@ss.bb.ba@izContainer@iZones[[1]]@angleGauge))
          stop('***>Angle gauge objects are not the same for Big BAF surfaces!')
        if(!identical(object@ss.ct.ba@izContainer@iZones[[1]]@angleGauge,
                      object@ss.ct.ba@izContainer@iZones[[1]]@angleGauge))
          stop('***>Angle gauge objects are not the same for count surfaces!')

        #check for the correct surface attributes...
        if(object@ss.bb.vol@estimate != .StemEnv$puaEstimates$volume ||
           object@ss.ct.vol@estimate != .StemEnv$puaEstimates$volume)
          stop('***>One of the volume surfaces is NOT a volume surface!')
        if(object@ss.bb.ba@estimate !=.StemEnv$puaEstimates$basalArea ||
           object@ss.ct.ba@estimate !=.StemEnv$puaEstimates$basalArea)
          stop('***>One of the basal area surfaces is NOT a basal area surface!')
          
        #make sure they all have the same number of stems at least...
        if(length(object@ss.bb.vol@izContainer@iZones) != 
               length(object@ss.bb.ba@izContainer@iZones) ||
           length(object@ss.bb.vol@izContainer@iZones) != 
               length(object@ss.ct.ba@izContainer@iZones) ||
           length(object@ss.bb.vol@izContainer@iZones) != 
               length(object@ss.ct.vol@izContainer@iZones))
          stop('***>All surfaces must have the same number of stems!')
          
        #they also must share the same tree population: very strict...
        strees.bb.vol = as(object@ss.bb.vol@izContainer, 'standingTrees')  #standingTrees container object
        strees.bb.ba = as(object@ss.bb.ba@izContainer, 'standingTrees')    #for BB basal area
        strees.ct.ba = as(object@ss.ct.ba@izContainer, 'standingTrees')    #for count basal area
        strees.ct.vol = as(object@ss.ct.vol@izContainer, 'standingTrees')  #for count volume
        if(!identical(strees.bb.vol, strees.bb.ba) || 
           !identical(strees.bb.vol, strees.ct.ba) ||
           !identical(strees.bb.vol, strees.ct.vol)
          )
          stop('***>"standingTree" objects must be exactly the same in all SS objects!')
        
                   
      return(TRUE)
    } #validity check
    
) #class ssBigBAF


