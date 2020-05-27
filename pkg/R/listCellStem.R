listCellStem = function(ss,
                        ...
                       )
{
#---------------------------------------------------------------------------
#
#   Perhaps a somewhat poor choice of name for this routine, it finds the
#   cell numbers that are within a stem's inclusion zone, and vice versa
#   (not within the actual stem!).
#
#   This routine was "lifted" from the /sampSurf/R/other file of the same
#   name. It has been slightly modified; I have also gotten rid of the 
#   option to make a SS object here. JHG, 1-May-2019.
#
#   This routine will take a sampling surface object and construct two
#   lists from the elements within...
#
#   1. A list by tree IDs of the grid cell numbers within each tree's
#      inclusion zone
#   2. A list of all of the grid cell numbers with elements the IDs for 
#      the stems whose inclusion zones overlap the individual 
#      cells--background cells are excluded since it's constructed from 
#      (1); essentially this is the "inverse" of (1)
#
#   Arguments...
#     ss = a "sampSurf" object
#     ... = gobbled
#
#   Returns...
#     a list invisibly with...
#       -- lists (1) & (2)
#
#   Note that this code can readily be used in sampSurf object construction
#   with a few changes, or after the fact from a SS object.
#
#Author...									Date: 12-Apr-2013 (original)
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   preliminaries: check for trees or logs, set things up...
#
    if(is(ss@izContainer, 'standingTreeIZs'))          #trees
      stems = as(ss@izContainer, 'standingTrees')
    else                                               #logs
      stems = as(ss@izContainer, 'downLogs')
    stemIDs = getID(stems)
    nstems = length(stemIDs)
    tract = ss@tract
 
#
#   get the list of cell numbers for each tree's inclusion zone...
#
    cellsInStems = vector('list', nstems)             #cell number by tree
    names(cellsInStems) = stemIDs
    for(i in 1:nstems ) {
      grid = izGrid(ss@izContainer@iZones[[i]], tract)@grid
      gvals = values(grid)
      idx = ifelse(is.na(gvals), FALSE, TRUE)         #index cells inside IZ, NAs are background
      gcells = seq_along(gvals)                       #all cell numbers
      gcells.iz = gcells[idx]                         #cell numbers within the IZ
      gxy = xyFromCell(grid, gcells.iz)               #xy values of the cell numbers
      egrid = extend(grid, tract)
      cells.iz = as.integer(cellFromXY(egrid, gxy))   #cell numbers in tract reference
      cellsInStems[[i]] = cells.iz                    #presumably (should be) in same order; could use id as key
    } #for

#
#   now "invert" it to list the stems in each cell--note that no background cells are
#   contained in this list, only cells within izones since it derives from the above...
#
    allCells = unlist(cellsInStems)
    uniqueCells = sort(unique(allCells))          #the unique set of cell numbers
    ncells = length(uniqueCells)
    stemsInCells = vector('list', ncells)
    names(stemsInCells) = as.character(uniqueCells)
    for(i in 1:nstems) {
      cdx = as.character(cellsInStems[[i]])
      for(j in cdx) 
#       stemsInCells[[j]] = c( stemsInCells[[j]], names(cellsInStems)[[i]] ) #no assigned standingTrees IDs
        stemsInCells[[j]] = c( stemsInCells[[j]], stemIDs[i])        #keeps stem ids from standingTrees too
    } #for

    return(invisible(list(cellsInStems = cellsInStems,
                          stemsInCells = stemsInCells
                         )
                    )
          )
}   #listCellStem
