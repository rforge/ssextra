getTreesBB = function(ssObj,
                      ...
                     )
{
#---------------------------------------------------------------------------
#
#   A simple utility to extract the trees from the sampling surface objects.
#
#   Note that it is general, and simply takes an object of class "sampSurf".
#
#   What makes it useful in BB sampling is the fact that it returns the
#   individual tree vbars too.
#
#   Arguments...
#     ssObj = a sampSurf object
#     ... = gobbled at present
#
#   Returns...
#     A list invisibly with...
#       -- a data frame of the trees in the population
#       -- a standingTrees object with all trees in the population
#       -- the number of trees common to each of the above
#
#   Updated for ssExtras and ssBigBAF, 26-Mar-2019, JHG.
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
#   check to be sure that what is passed is a sampSurf object...
#
    if(!is(ssObj, 'sampSurf'))
      stop('***>Object pass must be of class "sampSurf"!')
      
#
#   extract the standingTrees object and create a data frame augmented with 
#   volume and vbar...
#
    strees = as(ssObj@izContainer, 'standingTrees')               #standingTrees container object
    tree.ids = names(strees@trees)                                #generic tree ids for row names
    #tree.ids = trees$id                     #don't use the actual spatial tree ids for row names

#
#   now create a data frame from the above...
#
    nTrees = length(strees@trees)                                 #number of trees in the population
    trees = as(strees, 'data.frame')                              #cast to a data frame
    #add to the data frame...
    for(i in seq_len(nTrees)) {
      trees[i, 'id'] = getID(strees@trees[[i]])                   #tree spatial ids in this column
      trees[i, 'vol'] = strees@trees[[i]]@treeVol
      trees[i, 'ba'] = strees@trees[[i]]@ba
    }
    trees[, 'vbar'] = with(trees, vol/ba)
    rownames(trees) = tree.ids   
    
    return(invisible(list(trees = trees,
                          strees = strees,
                          nTrees = nTrees
                         )
                    )
          )
}   #getTreesBB


