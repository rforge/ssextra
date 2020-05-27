monteStatsBB = function(sdx = NA,
                        ssBB,
                        runQuiet = TRUE,
                        ...
                       )
{
#---------------------------------------------------------------------------
#
#   This routine will calculate the simple statistics (mean, var, var of mean)
#   on the individual sample draws from each of the four sampling surfaces
#   that can occur in big BAF sampling.
#
#   Note that the stats calculated here are over the grid cells/points in a
#   given sample draw. This is distinguished from the tree-based stats
#   that are calculated in monteTreeStatsBB().
#
#   Note: It may seem like a lot of overhead to pass the "ssBigBAF"
#         object. But my understanding of what happens is that only the
#         one (getValues) statement for each is evaluated as part of the
#         promise, and so, all that gets copied here are those values, not
#         the entire object. I tested it both ways (i.e., using getValues
#         in the calling routine and passing only theses values) with
#         microbenchmark and the difference was in the hundredths of a 
#         millisecond.
#
#   Arguments...
#     sdx = a vector of length m: the cell numbers in the sample
#     ssBB = an object of class "ssBigBAF"
#     runQuiet = TRUE: shhhh; FALSE: !shhhh
#     ... = gobbled at present
#
#   Returns...
#     A list invisibly with...
#       -- the data frame of samples
#       -- the data frame of stats for the count surfaces
#       -- the data frame of stats for the BB surfaces
#
#   Updates...
#     -- Updated for ssExtras and ssBigBAF, from surfStatsBB, 26-Mar-2019, JHG.
#     -- VBARS & Ratio vbar variance added 11-July-2019, JHG.
#
#
#Author...									Date: 12-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   check...
#
    if(!is(ssBB, 'ssBigBAF'))
      stop('***>Argument "ssBB" must be of class "ssBigBAF"!')

#
#   summarize for each surface...
#
    sdx = sort(sdx)
    cell.names = paste('c', sdx, sep = '.')
    
    n = length(sdx)                                         #sample size == number of sample points


#------------------------------------------
#   count (small BAF) surfaces (cs)...     
#------------------------------------------
#
    stats.ct = data.frame(matrix(NA_real_, nrow = 3, ncol = 2))
    colnames(stats.ct) = c('ba', 'vol')
    rownames(stats.ct) = c('mean', 'var', 'varm')

    #count ba surface..
    ba.cts = getValues(ssBB@ss.ct.ba@tract)[sdx]            #ba count: the sample of m cells/points
    #names(ba.cts) = cell.names
    stats.ct['mean', 'ba'] = mean(ba.cts)
    stats.ct['var', 'ba'] = var(ba.cts)
    stats.ct['varm', 'ba'] = stats.ct['var', 'ba']/n
    
    #count volume surface
    vol.cts = getValues(ssBB@ss.ct.vol@tract)[sdx]          #vol count: the sample of m cells/points
    #names(vol.cts) = cell.names
    stats.ct['mean', 'vol'] = mean(vol.cts)
    stats.ct['var', 'vol'] = var(vol.cts)
    stats.ct['varm', 'vol'] = stats.ct['var', 'vol']/n      #TGG notes (14); G&V (8.31)
    
    #the vbar stats are based on considering the ratio of vol to ba...
    mean.ba = stats.ct['mean', 'ba']
    mvbar.rat = stats.ct['mean', 'vol']/stats.ct['mean', 'ba']        #G&V (8.34) -- mean vbar
    S.r = sum( (vol.cts - mvbar.rat*ba.cts)^2 )/(n - 1)               #(8.43)
    r.varm = S.r/( (mean.ba * mean.ba) * n )                          #(8.42) var of the mean
    stats.ct['mean', 'vbar'] = mvbar.rat
    stats.ct['var', 'vbar'] = r.varm * n                              #not the variance of mean vbar
    stats.ct['varm', 'vbar'] = r.varm



#------------------------------
#   big BAF (bb) surfaces...   
#------------------------------
#
    stats.bb = data.frame(matrix(NA_real_, nrow = 3, ncol = 3))
    colnames(stats.bb) = c('ba', 'vol', 'vbar')
    rownames(stats.bb) = c('mean', 'var', 'varm')

    #bb vol surface...
    vol.bbs = getValues(ssBB@ss.bb.vol@tract)[sdx]          #vol bb: the sample of n cells/points
    #names(vol.bbs) = cell.names
    stats.bb['mean', 'vol'] = mean(vol.bbs)
    stats.bb['var', 'vol'] = var(vol.bbs)
    stats.bb['varm', 'vol'] = stats.bb['var', 'vol']/n      #TGG notes (14); G&V (8.31)
    
    #bb ba surface...
    ba.bbs = getValues(ssBB@ss.bb.ba@tract)[sdx]            #ba bb: the sample of n cells/points
    #names(ba.bbs) = cell.names
    stats.bb['mean', 'ba'] = mean(ba.bbs)
    stats.bb['var', 'ba'] = var(ba.bbs)
    stats.bb['varm', 'ba'] = stats.bb['var', 'ba']/n
    
    #the vbar stats are based on considering the ratio of vol to ba...
    mean.ba = stats.bb['mean', 'ba']
    mvbar.rat = stats.bb['mean', 'vol']/stats.bb['mean', 'ba']        #G&V (8.34) -- mean vbar
    S.r = sum( (vol.bbs - mvbar.rat*ba.bbs)^2 )/(n - 1)               #(8.43)
    r.varm = S.r/( (mean.ba * mean.ba) * n )                          #(8.42) var of the mean
    stats.bb['mean', 'vbar'] = mvbar.rat
    stats.bb['var', 'vbar'] = r.varm * n                              #not the variance of mean vbar
    stats.bb['varm', 'vbar'] = r.varm
  
#
#   now put all of the sample results into a data frame in case they are required;
#   ***Note that the cell values are totals for the tract size, since they come
#      directly from the cell values
#   the data frame is nrows = n cells (points)...
#
    samples = data.frame(sdx=sdx)                            #cell numbers
    #count-based...
    samples$ba.cts = ba.cts                                  #count basal area
    samples$vol.cts = vol.cts                                #count volume
    samples$idx.cts = ifelse(ba.cts > 0, TRUE, FALSE)        #index to count samples
    #big BAF-based...
    samples$ba.bbs = ba.bbs                                  #big BAF basal area
    samples$vol.bbs = vol.bbs                                #big BAF volume
    samples$idx.bbs = ifelse(vol.bbs > 0, TRUE, FALSE)       #index to BB samples
    #rownames(samples) = cell.names #these can be non-unique so don't use here
    
    return(invisible(list(samples = samples,
                          stats.ct = stats.ct,
                          stats.bb = stats.bb
                         )
                    )
          )
}   #monteStatsBB


