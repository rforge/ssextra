#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class for ssBigBAF class.
#
#   This uses the "sampSurf" class plotting routine, so please see the info
#   there on extensions; specifically the "sampSurf" plot routine uses the
#   "Tract" routine, which uses the "raster" routine...and so it goes...
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
#   method for class 'ssBigBAF'...
#
setMethod('plot',
          signature(x = 'ssBigBAF', y='missing'),
function(x,
         whichSS = slotNames(x)[grep('ss.', slotNames(x))],
         sampleCells = NULL,
         titles = NULL,
         namesAsTitles = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   This will plot one or more of the sampling surfaces in a "ssBigBAF" object.
#
#   Note that this will allow you to superimpose selected sample point locations
#   from say, a Monte Carlo draw, on top of the sampling surfaces. Simple pass
#   the cell numbers in the underlying raster in "sampleCells"...
#
#   Right now we have for whichSS...
#     1 = ss.bb.vol
#     2 = ss.bb.ba
#     3 = ss.ct.ba
#     4 = ss.ct.vol
#   this may be changed from numeric to the actual slot names above...
#
#   Here is an example where we are plotting the first Monte Carlo sample draw
#   from a set contained in m.p...
#
#   plot(z.p, whichSS=slotNames(z.p)[2:3], sampleCells=m.p@mc.samples$n.25$mc.1, 
#        point.col='red', namesAsTitles=TRUE)
#
#   note that "point.col" is passed to shwoCells.
#
#------------------------------------------------------------------------------
#
#   plot...
#
    ssNames = intersect(whichSS, slotNames(x)[1:4])
    np = length(ssNames)
    if(np == 0)
      stop('***>Choose surfaces by their slot names!')

    #decide on the arrangement based on the number of surfaces...
    apar = switch(as.character(np),
                    '1' = c(1,1),
                    '2' = c(1,2),
                    '3' = c(2,2),
                    c(2,2)
                   )
    oldpar = par(mfrow = apar)
    
    #this is the easy way out for titles: get it right or don't plot it...
    if(!is.null(titles) && length(titles) != np)
      titles = NULL
      
    #this usually does little good...
    suppressWarnings({                              #for non-plot arguments in ...
      
#
#   loop through each surface...
#
    for(i in seq_len(np)) {
      ss = paste0('x@', ssNames[i])
      plot.ss = parse(text=paste0('plot(', ss, ' ,...)'))
      eval(plot.ss)
      
      #plot the cells/points sampled if available...
      if(!is.null(sampleCells)) {
        SS = eval(parse(text = ss)) #extract the SS
        showCells(SS, cells = sampleCells, showSS = FALSE, ...)
      }
      
      #show title if there...
      if(!is.null(titles))
        title(titles[i])
      else if(namesAsTitles) #or default if desired
        title(ssNames[i])
    } #for plot loop

    })  #suppressWarnings 

    par(oldpar)
    return(invisible())

}    #plot for 'ssBigBAF'
) #setMethod


