initTract = function(extents = c(x = 150, y = 150),
                     cellSize = 1,                 
                     bufferWidth = 14,
                     units = 'metric',
                     ...
                    )
{
#---------------------------------------------------------------------------
#
#   Very simply, this sets up the Tract for the sampling surface simulations.
#
#   Arguments...
#     extents = the x,y extents relative from (0,0) for the tract
#     cellSize = the cell size in meters
#     bufferWidth = should be large enough to include the half-width of the
#                   largest inclusion zone for any sampling method tested
#     units = "English' or 'metric'
#     ... = other arguments that might be passed along to the Tract or
#           bufferedTract constructors
#
#   Returns...
#     a valid bufferedTract object invisibly
#
#   This routine has been handed down, most recently from the ComparingMethods
#   work.
#
#Author...									Date: 13-Mar-2017
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
    #require(sampSurf, quietly = TRUE)
    
    if(any(names(extents) != c('x','y')))
      stop('Extent names must be "x" and "y"!')
    
    tract = Tract(extents, cellSize=cellSize, units=units, ...)  #metric by default
    buffTr = bufferedTract(bufferWidth, tract, ...)              #buffer needs to include largest plot radius

    return(invisible( buffTr ))
} #initTract
 

