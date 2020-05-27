createBBNH = function(extents = c(x = 178, y = 178),
                      cellSize = 1,                 
                      bufferWidth = 18,  #to accommodate 3baf
                      units = 'metric',
                      baf.ct = 4,
                      baf.bb = 10,
                      startSeed = 355,
                      ...
                     )
{
#---------------------------------------------------------------------------
#
#   This little routine should be used to create a ssBigBAF object using
#   drawTreePop to generate a synthetic mixed northern hardwoods plot/stand.
#   The defaults will create the population used in the manuscript.
#
#   The number of trees drawn is determined from the basal area and the
#   Weibull parameters in drawTreePop--see it for more details.
#
#   A 142m x 142m tract with 18m buffer all around will create a tract that
#   is 3.17ha, or 7.83 acres; the actual inside-buffer plot is ~2ha in size
#   with a buffer area of 1.17ha. Thus, extents of 178x178 will produce the
#   tract with these dimensions. I.e., 142 + 2*18 = 178, to account for the
#   buffer on each side. This is all assuming a square grid cell size of 1m
#   on a side.
#
#   Arguments...
#     extents = the x,y extents relative from (0,0) for the tract
#     cellSize = the cell size in meters
#     bufferWidth = should be large enough to include the half-width of the
#                   largest inclusion zone for any sampling method tested
#     units = of measure
#     baf.ct = the basal area factor for the count sample
#     baf.bb = the basal area factor for the big BAF sample
#     startSeed = see jhgMisc::initRandomSeed
#     ... = passed on to other routines
#
#   Returns...
#     An object of class ssBigBAF invisibly
#
#Author...									Date: 10-June-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   create the tract first...
#
    #require(sampSurf, quietly = TRUE)

#
#   get the needed tract and tree population...
#
    nh = makePop(extents, cellSize, bufferWidth, units, startSeed = startSeed, ...)
    strees = nh$strees
    buffTr = nh$btr
    
#
#   create the standingTreeIZs inclusionZone container objects...
#
    cat('\nCreating standingtreeIZs...')
    ag.ct = angleGauge(baf.ct, units) 
    ct.strees = lapply(strees@trees, horizontalPointIZ, angleGauge = ag.ct)
    ct.izs = standingTreeIZs(ct.strees)
    
    ag.bb = angleGauge(baf.bb, units) 
    bb.strees = lapply(strees@trees, horizontalPointIZ, angleGauge = ag.bb)
    bb.izs = standingTreeIZs(bb.strees)

#
#   the individual sampling surfaces next...
#
    cat('\nCreating sampling surfaces...')
    cat('\n---Big BAF volume...')
    ss.bb.vol = sampSurf(bb.izs, buffTr, estimate = 'volume', ...)
    cat('\n---Big BAF basal area...')
    ss.bb.ba = sampSurf(bb.izs, buffTr, estimate = 'basalArea', ...)
    cat('\n---Count volume...')
    ss.ct.vol = sampSurf(ct.izs, buffTr, estimate = 'volume', ...)
    cat('\n---Count basal area...')
    ss.ct.ba = sampSurf(ct.izs, buffTr, estimate = 'basalArea', ...)

#
#   and finally the ssBigBAF object...
#
    cat('\nCreating ssBigBAF object...')    
    ssBB = ssBigBAF(ss.bb.vol = ss.bb.vol, ss.bb.ba = ss.bb.ba,            #big BAF
                    ss.ct.ba = ss.ct.ba, ss.ct.vol = ss.ct.vol,            #count BAF
                    ...)

    cat('\n***>Be careful to heed any warning about partial inclusion zones lying outside the tract!')
    cat('\n    This warning means you need to increase the buffer size for this population!')
    cat('\n')
    return(invisible(ssBB))


}   #createBBNH

