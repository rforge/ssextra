showCells = function(ss,
                     cells,
                     cell.col = 'salmon2',
                     alphaTrans = 0.2,
                     showSS = TRUE,
                     showPoints = TRUE,
                     point.col = 'gray25',
                     ...
                   )
{
#---------------------------------------------------------------------------
#
#   This will simply plot the cells given by the argument of the same name
#   onto a sampling surface to which it corresponds. It does this by creating
#   a new raster with extents of those in the bbox of the cells. It will
#   make any values in the raster NA that do not correspond to the cells
#   so that the non-NA values are the actual cell numbers and have nothing
#   to do with the values in the surface passed.
#
#   This is purely graphical, nothing more.
#
#   To display only the sample points, set showPoints=TRUE and alphaTrans=0.
#
#   Arguments...
#     ss = a sampSurf object
#     cells = a vector of cell numbers to be displayed
#     cell.col = the color for the "cells"
#     alphaTrans = ?transparentColorBase
#     showSS = TRUE: plot the sampSurf object; FALSE: assume we are adding
#              to an existing figure
#     showPoints = TRUE: show the sample points (cell centers); FALSE: don't
#     point.col = the color for the sample points
#     ... = gobbled at present
#
#   Returns...
#     The new raster object with "cells" only not NA
#
#   Here is an example were we simply draw a random sample out of a sampSurf
#   object that has 5625 cells...
#
#   showCells(ss.p10.vol, sample(1:5625, 100))
#
#Author...									Date: 8-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   should put some checks in here if this gets moved to sampSurf sometime...
#
    if(showSS)
      plot(ss)
      
#
#   in a new raster of the same size, assign only those cells in the "cells"
#   argument with non-zero (cell numbers) values, the rest are NA so that
#   they will not plot...
#
    cellNums = cellsFromExtent(ss@tract, extent(ss@tract))             #full set from ss
    rast = rasterFromCells(ss@tract, cellNums)                         #new raster, same extents
    r.vals = rep(NA_integer_, length(cellNums))                        #seed with NAs
    r.vals[cells] = cells                                              #cell numbers to show
    rast = setValues(rast, r.vals)                                     #replace
    
#
#   plot adds a legend, image will not...
#
    image(rast, add=TRUE, asp=1, col = transparentColorBase(cell.col, alphaTrans), ...)
#    plot(rast, add=TRUE, useRaster = FALSE,
#         col = transparentColorBase(cell.col, alphaTrans))

#
#   plot the sample points in the center of each cell if desired...
#
    if(showPoints) {
      xy = xyFromCell(rast, cells, spatial=TRUE)
      points(xy, pch=4, col = point.col, cex = 0.7, ...)
    }
    
    return(invisible(rast))
}   #showCells

