monteBootBB = function(samples,
                       B = 100,
                       alpha = 0.05,
                       ...
                      )
{
#---------------------------------------------------------------------------
#
#   This will perform bootstrap resampling as well as the jackknife on
#   a single MC sample drawn from the sampling surfaces. This routine
#   should be called from within the nested sample size within Monte Carlo
#   loops in monte(). The confidence intervals are computed simply
#   using the bcaboot::bcajack() routine.
#
#   It turns out that we can bootstrap (and jackknife) directly from the
#   samples component of the return from monteStatsBB() within the nested
#   loop in monte(). This contains the sample information at the different
#   samples$sdx cells/points in each of the four sampling surfaces for the
#   current MC draw on a given sample size, m. Some of the information in
#   samples is not need here, so we first trim that out below.
#
#   The bbMean() routine below does the estimation of the big BAF mean
#   for each bootstrap replication (also jackknife). One thing to note
#   here is that bcaboot::bcajack() does not seem to handle NAs at all.
#   Therefore, if we set tvbar.mean to NA_real_, as in monte() when
#   there are no trees sampled, we will get results of all NAs and NaNs
#   from the bootstrapping. Thus, here we set tvbar.mean = 0.0 in the
#   case of no sample trees.
#
#   Arguments...
#     samples = the samples data frame calculated in a single MC draw from
#               monteStatsBB() for a given sample size, m
#     B = the number of bootstrap samples
#     alpha = the two-tailed alpha level
#     ... = gobbled
#
#   Returns...
#     -- list of results from bcajack plus the jackknife estimate
#
#   Updated for ssExtras and ssBigBAF, from mcbootBB, 26-Mar-2019, JHG.
#
#Author...									Date: 19-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   we need only some of this, so make it more efficient...
#
    samples = samples[, c('ba.cts', 'ba.bbs', 'vol.bbs')]
    
#
#   this will calculate the big BAF volume estimate for each bootstrap sample...
#
    bbMean = function(x) {
              idx.bbs = ifelse(x$vol.bbs > 0, TRUE, FALSE)                     #index to BB samples
              mean.ba.cts = mean(x$ba.cts)
              tree.vbars = with(x, vol.bbs[idx.bbs]/ba.bbs[idx.bbs])
              tree.vbars = ifelse(is.nan(tree.vbars), NA_real_, tree.vbars)    #divide-by-zero check
              n.v = length(tree.vbars)
              tvbar.mean = ifelse(n.v > 0, mean(tree.vbars), 0.0)              #divide-by-zero check
              vol.bb = mean.ba.cts * tvbar.mean
              return(vol.bb)
    } #bbMean

#
#   limit the results by alpha, and bootstrap...
#
    alphas = c(alpha/2, 1-alpha/2)
    res = bcajack(samples, B, bbMean, alpha = alphas)

#
#   if desired, jackknife here--very wasteful since it is already computed by bcajack...
#
    jacksum = 0.0
    m = nrow(samples)
    for(i in seq_len(m)) {
      sx = samples[-i,]
      jacksum = jacksum + bbMean(sx)
    }
    jackmean = jacksum/m
    
    res$jackmean = jackmean
    
    return(res)
}   #monetBootBB
   
