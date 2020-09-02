#---------------------------------------------------------------------------
#
#   This augments an S4 method constructor for objects of class 
#   "monteBigBAF"...
#
#   Classes and subclasses in this file include...
#
#     1. monte  -- with object signature "monteBigBAF"
#
#   This was originally a stand-alone function that just operated on 
#   "ssBigBAF" and "monteBigBAF" objects as part of the PBDM paper
#   preparation. I have changed it a little to act like a constructor function
#   in the "monte" family.
#
#**>The "simplified" PBDM was added today. 1-Sept-2020, JHG.
#
#Author...									Date: 30-April-2020
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
# the generic for "monte" is already defined in sampSurf...
#
# thus, we must stick with the "object" argument below rather than change it
# to something more meaningful in this situation.


setMethod('monte',
          signature(object = 'monteBigBAF'),
#-----------------------------------------

function(object,
         ssBB = NA,
         runQuiet = FALSE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   This "monte" method will calculate Tom's PBDM variance approximation.
#   It uses simulation results from the main monte constructor and calculates
#   the new variance, covariance terms, correlation terms, and summary 
#   statistics (grand means, and the sampling variance of the MC means).
#
#   Note that "PBDM" stands for Point-Based Delta Method.
#
#   This method is based on the montePBDM.R version, which resides in the
#   R workspace code files for the second paper. However, this version has
#   been changed in that it calculates but does not save the "independent" 
#   PBDM (i.e., with no covariance terms) as this turned out to be very poor.
#
#   In the code the following are used as patterned after the TeX code in
#   the manuscript...
#      Bc = Basal area from the count points
#      Bv = Basal area from the big BAF (i.e., volume) points
#      Vv = Volume from the big BAF points
#
#   Arguments...
#     object = the corresponding "monteBigBAF" object (ssBB.mc)
#     ssBB = an "ssBigBAF" object used in the original simulations
#     runQuiet = TRUE: no feedback; FALSE: a short report on the results
#     ... =  gobbled for now
#
#   Returns...
#     a "monteBigBAF" object invisibly that has been augmented with PBDM stats
#
#Author...									Date: 15-Nov-2019 prototype
#	Jeffrey H. Gove                               29-April-2020 current
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   make sure we have routines to work with...
#
    #require(sampSurf, quietly = TRUE)
    
    ssBB.mc = object
    if(!is(ssBB, 'ssBigBAF'))
      stop('ssBB object not of class "ssBigBAF"')

#
#   get the MC sample and sample size info...
#
    mcSamples = ssBB.mc@mcSamples                     #number of MC samples
    mc.names = names(ssBB.mc@mc.samples[[1]])         #they will be the same for any sample size, n
    n.mcs = length(mc.names)

    n = ssBB.mc@n
    n.names = names(n)
    nn = length(n)

#
#   t-values for alpha level desired...
#
    alpha = ssBB.mc@alpha
    alpha2 = alpha/2
    t.values = qt(1-alpha2, n - 1)
    names(t.values) = n.names
    
    popVolume = ssBB@ss.bb.vol@surfStats$popTotal  #this will be the same for ss.ct.vol
    
#
#   set up the storage lists (of data frames) for the results: patterned after
#   the ssExtra::monte constructor for "ssBigBAF" objects...
#
    rvec = data.frame(matrix(NA_real_, nrow = n.mcs, ncol = nn))
    rownames(rvec) = mc.names
    colnames(rvec) = n.names
    
    #variances of the means...
    varMeans = list(pbDelta.ind = rvec,     #assuming independence (just to have this recorded)
                    pbDelta = rvec,         #full PBDM with covariances
                    pbdelta = rvec          #simplified PBDM
                   )
                   
    #standard error of the means...
    stErrs = varMeans
    
    #associated lower and upper confidence intervals limits...
    lowerCIs = varMeans
    upperCIs = varMeans
                
    #covariances of the means...
    covs = list(BcVv.cov = rvec,         #Count BA (Bc) vs. BB vol (Vv)
                BcBv.cov = rvec,         #Count BA (Bc) vs. BB BA (Bv)
                VvBv.cov = rvec          #BB vol (Vv) vs. BB BA (Bv)
               )
               
    #correlations...
    cors = list(BcVv.cor = rvec,
                BcBv.cor = rvec,
                VvBv.cor = rvec
               )
               
    #whether the CIs caught the population mean (TRUE) or no (FALSE)...
    lvec = rvec
    lvec[] = FALSE
    caught = list(pbDelta.ind = lvec,        #all as in varMeans, but with logical seed
                  pbDelta = lvec,
                  pbdelta = lvec
                 ) 

               
#
#-----------------------------------------------------------------------------
#   Loop through each Monte Carlo sample (i) within sample size (j)...
#
    if(!runQuiet) cat('\nn=')
    for(j in seq_len(nn)) {                               #sample size first
      n.j = n[j]
      n.lab = n.names[j]
      if(!runQuiet)
        cat(n.j,', ',sep='')
      for(i in seq_len(n.mcs)) {
        #calculate the point-based stats on the current sample of cells (points)...
        mc.lab = mc.names[i]
        ms = monteStatsBB(ssBB.mc@mc.samples[[n.lab]][[mc.lab]], ssBB)
 
        #BB stats...
        stats.bb = ms$stats.bb
        Vv = stats.bb['mean', 'vol']              #volume from the BB points
        Bv = stats.bb['mean', 'ba']               #basal area from the BB points
        Vv.varm = stats.bb['varm', 'vol']
        Bv.varm = stats.bb['varm', 'ba']
    
        #count stats...
        stats.ct = ms$stats.ct
        Bc = stats.ct['mean', 'ba']               #basal area from the count points
        Bc.varm = stats.ct['varm', 'ba']

#
#       the derivative multipliers and 'independent' component of the full variance...
#
        Vv.Bv = Vv/Bv
        Bc.Bv = Bc/Bv
        VvBc.Bv = Vv*Bc/(Bv^2)
    
        pbDelta.ind = Vv.Bv^2 * Bc.varm + Bc.Bv^2 * Vv.varm + VvBc.Bv^2 * Bv.varm
        varMeans$pbDelta.ind[mc.lab, n.lab] =  pbDelta.ind
        stErrs$pbDelta.ind[mc.lab, n.lab] = sqrt(pbDelta.ind)
        
#
#       calculate the covariances (of the means) and the full approximate variance...
#
        samples = ms$samples
        BcVv.cov = with(samples, cov(ba.cts, vol.bbs)/n.j)
        BcBv.cov = with(samples, cov(ba.cts, ba.bbs)/n.j)
        VvBv.cov = with(samples, cov(vol.bbs, ba.bbs)/n.j)
        pbDelta = pbDelta.ind + 2*Vv.Bv*Bc.Bv * BcVv.cov - 
                                2*VvBc.Bv*Vv.Bv * BcBv.cov -
                                2*VvBc.Bv*Bc.Bv * VvBv.cov
        varMeans$pbDelta[mc.lab, n.lab] = pbDelta
        stErrs$pbDelta[mc.lab, n.lab] = sqrt(pbDelta)
                                      
        covs$BcVv.cov[mc.lab, n.lab] = BcVv.cov
        covs$BcBv.cov[mc.lab, n.lab] = BcBv.cov
        covs$VvBv.cov[mc.lab, n.lab] = VvBv.cov
        
#
#       the Expectation approximation in the paper based on P&H simplification...
#
        BcVv = Bc*Vv
        BcBc = Bc*Bc
        VvVv = Vv*Vv
        pbdelta = VvVv*(Bc.varm/BcBc + Vv.varm/VvVv + Bv.varm/BcBc +
                            2*( BcVv.cov/BcVv - BcBv.cov/BcBc -
                                VvBv.cov/BcVv ) #cov terms
                           )
        varMeans$pbdelta[mc.lab, n.lab] = pbdelta
        stErrs$pbdelta[mc.lab, n.lab] = sqrt(pbdelta)

#
#       correlations...
#
        BcVv.cor = BcVv.cov/(sqrt(Bc.varm)*sqrt(Vv.varm))
        BcBv.cor = BcBv.cov/(sqrt(Bc.varm)*sqrt(Bv.varm))
        VvBv.cor = VvBv.cov/(sqrt(Vv.varm)*sqrt(Bv.varm))

        cors$BcVv.cor[mc.lab, n.lab] = BcVv.cor
        cors$BcBv.cor[mc.lab, n.lab] = BcBv.cor
        cors$VvBv.cor[mc.lab, n.lab] = VvBv.cor
        
#
#       confidence intervals...
#
        x = .ciCaught(ssBB.mc@means$vol.bb[i,j], stErrs$pbDelta.ind[mc.lab, n.lab], 
                     t.values[j], popVolume)
        caught$pbDelta.ind[i,j] = x$catch
        lowerCIs$pbDelta.ind[i,j] = x$lci
        upperCIs$pbDelta.ind[i,j] = x$uci

        x = .ciCaught(ssBB.mc@means$vol.bb[i,j], stErrs$pbDelta[mc.lab, n.lab], 
                     t.values[j], popVolume)
        caught$pbDelta[i,j] = x$catch
        lowerCIs$pbDelta[i,j] = x$lci
        upperCIs$pbDelta[i,j] = x$uci

        x = .ciCaught(ssBB.mc@means$vol.bb[i,j], stErrs$pbdelta[mc.lab, n.lab], 
                     t.values[j], popVolume)
        caught$pbdelta[i,j] = x$catch
        lowerCIs$pbdelta[i,j] = x$lci
        upperCIs$pbdelta[i,j] = x$uci

      } #mcSampleNos for loop
    } #sample size loop


#
#   save the above results to the monteBigBAF object (not the *.ind results)...
#
    ssBB.mc@varMeans = c(ssBB.mc@varMeans, varMeans['pbDelta'], varMeans['pbdelta'])
    ssBB.mc@stErrs = c(ssBB.mc@stErrs, stErrs['pbDelta'], stErrs['pbdelta'])
    ssBB.mc@corrs = c(ssBB.mc@corrs, cors['BcVv.cor'])
    ssBB.mc@corrs = c(ssBB.mc@corrs, cors['BcBv.cor'])
    ssBB.mc@corrs = c(ssBB.mc@corrs, cors['VvBv.cor'])
    ssBB.mc@covs = c(ssBB.mc@covs, covs['BcVv.cov'])
    ssBB.mc@covs = c(ssBB.mc@covs, covs['BcBv.cov'])
    ssBB.mc@covs = c(ssBB.mc@covs, covs['VvBv.cov'])
    ssBB.mc@caught = c(ssBB.mc@caught, caught['pbDelta'], caught['pbdelta'])
    ssBB.mc@lowerCIs = c(ssBB.mc@lowerCIs, lowerCIs['pbDelta'], lowerCIs['pbdelta'])
    ssBB.mc@upperCIs = c(ssBB.mc@upperCIs, upperCIs['pbDelta'], upperCIs['pbdelta'])


#
#   average summary statistics--grand means over each sample for each n;
#   note that it also can contain jackknife & bootstrap results...
#
    gm.all = list(gm.varMeans = .grandMeans(varMeans, n.names),
                  gm.stErrs = .grandMeans(stErrs, n.names),
                  gm.covs = .grandMeans(covs, n.names),
                  gm.cors = .grandMeans(cors, n.names),
                  gm.caught = .grandMeans(caught, n.names) * 100.0,
                  gm.lowerCIs = .grandMeans(lowerCIs, n.names),
                  gm.upperCIs = .grandMeans(upperCIs, n.names)
                 )
    
#
#   save to the monteBigBAF object...
#
    ssBB.mc@gm.all$gm.varMeans = rbind(ssBB.mc@gm.all$gm.varMeans, 
                                       gm.all$gm.varMeans['pbDelta', , drop=FALSE],
                                       gm.all$gm.varMeans['pbdelta', , drop=FALSE]
                                     )
    ssBB.mc@gm.all$gm.stErrs = rbind(ssBB.mc@gm.all$gm.stErrs, 
                                     gm.all$gm.stErrs['pbDelta', , drop=FALSE],
                                     gm.all$gm.stErrs['pbdelta', , drop=FALSE]
                                    )
    ssBB.mc@gm.all$gm.corrs = rbind(ssBB.mc@gm.all$gm.corrs, gm.all$gm.cors)
    ssBB.mc@gm.all$gm.caught = rbind(ssBB.mc@gm.all$gm.caught, 
                                     gm.all$gm.caught['pbDelta', , drop=FALSE],
                                     gm.all$gm.caught['pbdelta', , drop=FALSE]
                                    )
    ssBB.mc@gm.all$gm.lowerCIs = rbind(ssBB.mc@gm.all$gm.lowerCIs, 
                                       gm.all$gm.lowerCIs['pbDelta', , drop=FALSE],
                                       gm.all$gm.lowerCIs['pbdelta', , drop=FALSE]
                                      )
    ssBB.mc@gm.all$gm.upperCIs = rbind(ssBB.mc@gm.all$gm.upperCIs, 
                                       gm.all$gm.upperCIs['pbDelta', , drop=FALSE],
                                       gm.all$gm.upperCIs['pbdelta', , drop=FALSE]
                                      )
     
    
#
#-----------------------------------------------------------------------------
#   the approximate sampling variances of the means from the MC replications;
#   computed for each of the four sampSurf objects...
#
#    sm.all = .samplingVarMeans(ssBB.mc)       #This is now calculated in the main monte<<******

    if(!runQuiet) 
      cat('\n')

#
#   send it back...
#
    return(invisible(ssBB.mc))
    
}   #montePBDM
)   #setMethod

