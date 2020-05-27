makePop = function(extents = c(x = 150, y = 150),
                   cellSize = 1,                 
                   bufferWidth = 14,
                   units = 'metric',
                   description = "Tract and northern hardwoods population",
                   ...
                  )
{
#---------------------------------------------------------------------------
#
#   Very simply, this makes the tract and tree population that are used
#   in the simulations, just take the defaults.
#
#   Note that the above defaults will create a buffered tract with total
#   area of 2.25ha and an area internal to the buffer of ~1.5ha (1.4884).
#   These are the defaults used for the manuscript.
#
#   This...
#     -- makes the buffered tract using initTract()
#     -- calls drawTreePop() to create the tree population on the tract
#     -- creates a "standingTrees" version of the tree pop data frame
#
#   Arguments...
#     extents = the x,y extents relative from (0,0) for the tract
#     cellSize = the cell size in meters
#     bufferWidth = should be large enough to include the half-width of the
#                   largest inclusion zone for any sampling method tested
#     ... = passed on to drawTreePop()
#
#   Returns...
#     a list invisibly with...
#       -- the buffered tract object
#       -- the standingTrees object
#       -- corresponding trees data frame
#
#Author...									Date: 29-Mar-2017
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   tract first...
#
    #require(sampSurf, quietly = TRUE)
 
    cat('\nCreating tract...')   
    btr = initTract(extents = extents, cellSize = cellSize, bufferWidth = bufferWidth,
                   units = units, description = description) #no ... for replication purposes
    
#
#   now the tree population...
#
    cat('\nCreating tree population...')
    trees = drawTreePop(btr, showPlot = FALSE, ...)
    cat('\nCreating standingTrees object...')
    strees = standingTrees(trees)

#
#   add ba, vbars and spIDs to the data frame...
#
    cat('\nAdding to data frame...')
    tree.ids = names(strees@trees)                                #generic tree ids for row names
    nTrees = length(strees@trees)                                 #number of trees in the population
    for(i in seq_len(nTrees)) {
      trees[i, 'id'] = getID(strees@trees[[i]])                   #tree spatial ids in this column
      trees[i, 'vol'] = strees@trees[[i]]@treeVol
      trees[i, 'ba'] = strees@trees[[i]]@ba
    }
    trees[, 'vbar'] = with(trees, vol/ba)
    rownames(trees) = tree.ids   

    
    cat('\n')
    return(invisible(list(btr = btr,
                          strees = strees,
                          trees = trees
                         )
                    )
          )
    
    
}   #makePop


