ratioVariance = function(ssBB,
                         ssBB.mc,
                         n = 25,
                         mcSample = 1,
                         stype = 'bigBAF',
                         runQuiet = FALSE,
                         ...
                        )
{
#---------------------------------------------------------------------------
#
#   This routine will calculate the ratio variance for the volume/ba ratio
#   in (8.34) of G&V. The variance is calculated via (8.42) & (8.43); both
#   being found on p. 258-259.
#
#   Note that all of what is here is in the monteBigBAF constructor, but not
#   via this routine, which is simply made available to calculate the
#   variances for individual MC draws if desired, outside the framework
#   of the constructor.
#
#   Also note that when this was written, the ratio variance estimation was
#   not part of the monteStatsBB() routine. Therefore, the calculations were
#   performed here: this routine was a prototype and the code below was
#   ported to monteStatsBB.
#
#   Arguments...
#     ssBB = an object of class 'ssBigBAF'
#     ssBB.mc = an object of class 'monteBigBAF'
#     n = the sample size == number of sample points/cells
#     mcSample = the MC sample number to use given n
#     stype = either 'bigBAF' or 'count' -- see monteTreeStatsBB()
#     runQuiet = TRUE: ssshhh; FALSE: hmmm
#     ... = gobbled
#
#   Returns...
#     a list invisibly with...
#       -- the results of monteStatsBB() for n & mcSample
#       -- the results of monteTreeStatsBB() for n & mcSample
#       -- the sample size
#       -- the MC sample number
#       -- a data frame with the variance of the mean and standard
#          errors for the ratio and normal theory estimators
#       -- a vector of the means
#
#Author...									Date: 9-July-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   checks & setup...
#
    if(!is(ssBB, 'ssBigBAF'))
      stop('ssBB argument must be of class "ssBigBAF"!')
    if(!is(ssBB.mc, 'monteBigBAF'))
      stop('ssBB.mc argument must be of class "monteBigBAF"!')
    if(!n %in% ssBB.mc@n)
      stop('illegal sample size, n, not found in MC sims!')
    if(mcSample < 1 || mcSample > ssBB.mc@mcSamples)
      stop(paste('mcSample must be in the range of [',1, ',', ssBB.mc@mcSamples, ']', sep=''))
      
    n.ss = paste('n', n, sep='.')                          #e.g., "n.25"
    mc.num = paste('mc', mcSample, sep='.')                #e.g., "mc.1"
    mc.sample = ssBB.mc@mc.samples[[n.ss]][,mc.num]        #a list and a dataframe

#
#   get the stats for plot-level and tree-level...
#
    msBB.mc = monteStatsBB(mc.sample, ssBB)
    mtsBB.mc = monteTreeStatsBB(mc.sample, ssBB, stype)
     

#
#   collect the required quantities for the ratio variance...
#
    if(stype == 'bigBAF') {
      mvbar.rat = with(msBB.mc$stats.bb['mean',], vol/ba)                   #G&V (8.34)
      mean.vbar = mtsBB.mc$treeStats['mean', 'vbar']                        #mean vbar from trees
      vol.pws = msBB.mc$samples$vol.bbs                                     #point-wise volume
      ba.pws = msBB.mc$samples$ba.bbs                                       #point-wise ba
      npts = length(ba.pws)                                                 #no of sample points
      mean.ba = msBB.mc$stats.bb['mean', 'ba']
    }
    else {                                                                  #count sample case...
      mvbar.rat = with(msBB.mc$stats.ct['mean',], vol/ba)                   #G&V (8.34)
      mean.vbar = mtsBB.mc$treeStats['mean', 'vbar']                        #mean vbar from trees
      vol.pws = msBB.mc$samples$vol.cts                                     #point-wise volume
      ba.pws = msBB.mc$samples$ba.cts                                       #point-wise ba
      npts = length(ba.pws)                                                 #no of sample points
      mean.ba = msBB.mc$stats.ct['mean', 'ba']
    }
 
#
#   and the ratio variance estimates...
#
    S.r = sum( (vol.pws - mean.vbar*ba.pws)^2 )/( npts*(npts - 1) )       #(8.43)/n
    r.varm = S.r/(mean.ba * mean.ba)                                      #(8.42) var of the mean
    r.se = sqrt(r.varm)
    
#
#   a check with the monte results...
#
    if(stype == 'bigBAF')
      R.varm = msBB.mc$stats.bb['varm', 'vbar']
    else
      R.varm = msBB.mc$stats.ct['varm', 'vbar']
    
    nt.varm = mtsBB.mc$treeStats['varm', 'vbar']                          #normal theory var of mean

#
#   a little feedback if desired...
#
    if(!runQuiet) {
      cat('\nSample size n =', n)
      cat('\nMonte Carlo replicate =', mcSample)
      cat('\nRatio variance of the mean =', r.varm, 'calculated here')
      cat('\nRatio variance of the mean =', R.varm, 'from monteStatsBB')
      cat('\nNormal theory varm =', nt.varm)
      cat('\n')
    }

#
#   just compare what was calculated here with the monteStatsBB results...
#
    aet = all.equal(R.varm, r.varm) 
    if(is.character(aet))
      cat('\nRatio variance comparison:', aet)

#
#   wrap it up & return...
#
    df = data.frame(matrix(NA_real_, nrow = 2, ncol = 2))
    rownames(df) = c('varm', 'se')
    names(df) = c('ratio', 'normTheory')
    
    df['varm', 'ratio'] = r.varm
    df['varm', 'normTheory'] = nt.varm
    df['se', 'ratio'] = sqrt(r.varm)
    df['se', 'normTheory'] = sqrt(nt.varm)

    return(invisible(list(msBB.mc = msBB.mc,
                          mtsBB.mc = mtsBB.mc,
                          n = n,
                          mcSample = mcSample,
                          df = df,
                          means = c(mean.ba = mean.ba, mean.vbar = mean.vbar, mvbar.rat = mvbar.rat) 
                         )
                    )
          )

}   #ratioVariance

