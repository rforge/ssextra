monteTreeStatsBB = function(sdx,
                            ssBB,
                            stype = c('bigBAF', 'count'),
                            ...
                           )
{
#---------------------------------------------------------------------------
#
#   This routine performs a similar function to wtreeVbars, but where the
#   wtreeVbars uses the csl@cellsInStems list to get its results, this
#   routine uses the csl@stemsInCells list. The approach used here appears
#   to be significantly faster as demonstrated below and thus supercedes 
#   wtreeVbars.
#
#   Note: wtreeVbars is not included in this package, it is found in...
#         ./other/wtreeVbars.R in this directory. It was used in some of
#         the prototype code found in: /home/jhgove/ForestInventory/ #+
#         HorizontalPointSampling/DoubleSampling/BigBAF/Variance/Rwork/R
# 
#   In addition, either the big BAF-based or the Count-based summaries
#   can be selected w/ the "stype" argument. Thus, if only one is needed,
#   this cuts down on computer time. If both are needed, running it 
#   twice is still quicker than running wtreeVbars.
#
#   The important point to remember w/r to this routine is that all the 
#   results are tree-level summaries, nothing is expanded. This is what
#   distinguishes it from monteStatsBB(), which returns stats based on
#   cell/point-wise expanded totals, which is why vbar is a major component
#   in the current routine, while it is missing in monteStatsBB().
#
#   Note especially then that here ba is the actual tree basal area for
#   each tree (in df.trees) or the sum by point/cell (in df.cells).
#   Thus, to match the output in terms of ba from monteStatsBB()
#   in samples$ba.bbs (for BB points) or samples$ba.cts (for count points)
#   simply multiply df.cells$no.trees*F*areaAdjust, where
#       -- F is the basal area factor
#       -- areaAdjust = area(tract)/unitArea (see the sampSurf constructor)
#
#   Here are some comparison stats (the arguments below were a little
#   different than in the current version)...
#
#     microbenchmark({monteTreeStatsBB(mc410.p6@mc.samples$n.25$mc.1, 
#                     ss410.p@csl.bb, df$trees)})
#            min        lq      mean    median        uq       max neval
#      10.632416 10.831032 11.252953 11.006056  11.30816 15.249227   100
#
#     microbenchmark({wtreeVbars(mc410.p6@mc.samples$n.25$mc.1, 
#                     ss410.p@csl.bb@cellsInStems,  
#                     ss410.p@csl.ct@cellsInStems, df$trees)})
#             min        lq      mean    median        uq       max neval
#       28.558108 29.656691 31.192882 30.264909 32.668487 43.824179   100
#
#   and just for comparison between bb and count...
#
#     microbenchmark({monteTreeStatsBB(mc410.p6@mc.samples$n.25$mc.1, 
#                     ss410.p@csl.ct, df$trees)})
#             min        lq      mean    median        uq       max neval
#       10.759261 10.831825 11.164804 10.887148 11.025626 14.537476   100
#
#   Thus, this routine is approximately three times faster than wtreeVbars()
#   and this routine is also calculating basal area & stats. The count
#   results are comparable, even though the number of sampled trees is larger.
#
#   Arguments...
#     sdx = the sample cell numbers; sdx = sample point index originally--
#           note that it will be renamed to cellNum in the data frame
#           exported here; this will typically be a sample of cells from
#           monteBBConstructor()
#     ssBB == an object of class"ssBigBAF"
#     stype = "bigBAF" for tree stats on the big BAF surfaces; "count" for
#             the count surfaces
#     ... = gobbled
#
#   Returns...
#     A list invisibly with...
#       -- a data frame based on the trees measured; but it does include
#          records for background sample points, so tree-based stats
#          whould be run using the idx column to pick out sampled cells;
#          this concatenates the list of data frames below...
#       -- a list of data frames corresponding to the above, where the
#          individual data frames are NULL for packground points or
#          the trees measured in that cell/point; this will have length
#          length(sdx), or the sample size
#       -- a cell/point-based data frame where the trees have been 
#          aggregated (summed) by cell for totals at each cell/point
#       -- a data frame with tree-based summary statistics
#       -- a data frame with cell/point-based summary statistics
#       -- the number of sampled vbar trees
#
#**>Note above that the ssBB@csl.* slot contains the stemsInCells list for
#   either the big BAF or Count surfaces. Now, as constructed, these lists
#   contain ONLY the non-background cells. Thus, in the code below, the
#   creation of the initial sampled list will make NULL list entries for
#   background cells, since they don't exist in the stemsInCells list. This
#   is convenient, and allows us to preserve this record of non-tally points
#   in the final full data frames.
#
#**>Note: when sampling with replacement, cell numbers can occur in the
#   stemsInCells list multiple times. Serendiptously, the list structure
#   aids us in keeping these unique sample points separated via the
#   structure of having one list element for each cell/point. A unique
#   but meaningless "replicate/replacement" cell index is added to the trees
#   data frame below to allow us to aggregate by this index first, then cell
#   number within the index for the cell/point-wise summary combining
#   all trees by points. Were aggregation only over cell number, it would
#   not keep replicated visits to the same cell (replacement) separate, they
#   would be aggregated into one. This mistake would yield a cells df with
#   too few sample points. This rep.cdx (cell index) is also found in the
#   df.cells aggregated data frame, but because of aggregation by point
#   there will be no repeats in the rep.cdx in this data frame, where there
#   can be repeats in the df.trees data frame. Make sense?
#
#Author...									Date: 6-May-2019+
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   extract the data frame of trees with vbars & see if we are going to
#   process BB or Count samples; since the trees are all the same in
#   each of the four sampling surfaces in an ssBigBAF object, it does not 
#   matter which one is used to get the getTreesBB()$trees below...
#
    trees = getTreesBB(ssBB@ss.bb.vol)$trees    #same trees on each surface
    stype = match.arg(stype)
    if(stype == 'bigBAF')                       #are we doing BB or count?
      csl = ssBB@csl.bb
    else
      csl = ssBB@csl.ct
#
#   get the subset by plot number; see the note above, all background cells
#   will automatically receive a list entry with the cell number as a name
#   but with NULL contents...
#
    sdx = sort(sdx)                       #sort to be identical to monteStatsBB()
    cdx = as.character(sdx)               #cell/point numbers, including repeats if w/ replacement
    stemsInCells = csl@stemsInCells       #full surface list of stem numbers at each cell/point
    sic = stemsInCells[cdx]               #we want only those cells in the sample sdx
    names(sic) = cdx                      #some will most assuredly be NA (background)
    n = length(sic)                       #the sample size

#
#   the following is a unique index that is required to keep sample points with the sample cell
#   number from being aggregated together in the df.cells data frame below;
#   this is necessary for sampling with replacement where, for example, two background cells
#   can be sampled twice as an easy illustration; also two cells inside tree inclusion zones
#   could be sampled independently with the same tree list, and these must not be aggregated
#   into one record; this independent index added to the df.trees below will assure that
#   the multiply sampled replacement point/cells will remain separate...
#
    rdx = seq_along(cdx)                  #quite simple: a unique id for each sample point

#
#   a little function to determine whether a list element is NULL or not and returns
#   a newly created data frame including the id, vbar and ba for the matching trees
#   from the "trees" data frame...
#
    pick = function(x) {
             if(is.null(x)) 
               return(x)
               #return(c(NA_real_, 0.0, 0.0))
             else {
               idx = match(unlist(x), trees$id)                                   
               return(trees[idx, c('id', 'vbar', 'ba')])
             }
           } #pick
    
#
#   create a list of data frames for each cell/point with one or more trees; or NULL for
#   no sampled trees on the point...
#
    list.df = lapply(sic, pick)

#
#   the null prototype data frame for use in constructing the full df below...
#
    df.null = data.frame(rep.cdx = NA_integer_,     #replicated cell index--sampling with replacement
                         cellNum = NA_integer_,     #cell (point) number in the surface
                         tree.id = '',              #the relative tree number (rownames)
                         id = '',                   #the tree's original spatial id: spID
                         vbar = 0.0,                #hmmmm
                         ba = 0.0,                  #hmmmmmmmm
                         idx = FALSE                #FALSE: background; TRUE: sampled
                        )   

#
#   loop through all cells/points and cbind the cell data frames from list.df into one; note that
#   the result will be a data frame with individual records for every tree sampled on each
#   cell/point; thus, if a tree is sampled on more than one point, it will show up multiple times
#   in the data frame...
#
    df.trees = NULL                      #dummy to begin for easy concatenation
    for(i in seq_len(n)) {
      if(is.null(list.df[[i]])) {        #background: dummy record with no trees for the point
        df.null$rep.cdx = rdx[i]
        df.null$cellNum = sdx[i]
        df.null$id = NA_character_
        df.null$tree.id = NA_character_
        df.trees = rbind(df.trees, df.null)
      }
      else { 
        rep.no = rep(rdx[i], nrow(list.df[[i]]))
        cell.no = rep(sdx[i], nrow(list.df[[i]]))
        zz = cbind(rep.cdx = rep.no, cellNum = cell.no, 
                   tree.id = rownames(list.df[[i]]), list.df[[i]], 
                   idx = TRUE
                  )
        df.trees = rbind(df.trees, zz)
      }
    } #for
    rownames(df.trees) = seq(1, nrow(df.trees))

#
#   this aggregates the vbars and ba over all trees on each cell/point for point totals;
#   note that this will include zero-tally (background) cells/points; also note that
#   aggregating idx give us the number of trees...
#
#   importantly, this aggregates by rep.cdx then cell number within rep.cdx so
#   that replicate ("replacement") cells will be kept separated in the aggregation...
#
    df.cells = aggregate(df.trees[,c('vbar','ba', 'idx')], 
                         list(rep.cdx = df.trees$rep.cdx, cellNum = df.trees$cellNum), 
                         sum
                        )
    names(df.cells)[5] = 'no.trees'   #rename idx so it won't get overwritten below!
    df.cells$idx = ifelse(df.cells$vbar > 0.0, TRUE, FALSE)     #just in case we need it
    n.v = sum(df.cells$no.trees)                        #number of trees measured for vbars

#
#   tree-based stats: based on n.v trees, ignore background...
#
    treeStats = data.frame(matrix(NA_real_, nrow = 4, ncol = 2))
    names(treeStats) = c('vbar', 'ba')
    rownames(treeStats) = c('mean', 'var', 'varm', 'corr')
    treeStats['mean', 'vbar'] = with(df.trees, mean(vbar[idx]))  
    treeStats['var', 'vbar'] = with(df.trees, var(vbar[idx]))
    treeStats['varm', 'vbar'] = treeStats['var', 'vbar']/n.v
    treeStats['mean', 'ba'] = with(df.trees, mean(ba[idx]))  
    treeStats['var', 'ba'] = with(df.trees, var(ba[idx]))
    treeStats['varm', 'ba'] = treeStats['var', 'ba']/n.v
    treeStats['corr', ] = with(df.trees, cor(vbar[idx], ba[idx]))          #tree vbar-ba by tree

#
#   cell/point-based summary statistics; note that the cell/point-based
#   stats DO include the background/zero cells since they are part of 
#   the point-wise sample of n cells...
#
    cellStats = data.frame(matrix(NA_real_, nrow = 4, ncol = 2))
    names(cellStats) = c('vbar', 'ba')
    rownames(cellStats) = c('mean', 'var', 'varm', 'corr')
    cellStats['mean', 'vbar'] = with(df.cells, mean(vbar))  
    cellStats['var', 'vbar'] = with(df.cells, var(vbar))
    cellStats['varm', 'vbar'] = cellStats['var', 'vbar']/n
    cellStats['mean', 'ba'] = with(df.cells, mean(ba))  
    cellStats['var', 'ba'] = with(df.cells, var(ba))
    cellStats['varm', 'ba'] = cellStats['var', 'ba']/n
    cellStats['corr', ] = with(df.cells, cor(vbar, ba))             #aggregate tree vbar-ba


#
#   send it all back...
#
    return(invisible(list(df.trees = df.trees,
                          list.df = list.df,
                          df.cells = df.cells,
                          treeStats = treeStats,
                          cellStats = cellStats,
                          n.v = n.v
                         )
                    )
          )
}   #monteTreeStatsBB

