#---------------------------------------------------------------------------
#
#   This file holds the S4 method constructor for objects of class 
#   "monteBigBAF"...
#
#   Classes and subclasses in this file include...
#
#     1. monte  -- with object signature "ssBigBAF"
#
#Author...									Date: 20-Mar-2019
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
          signature(object = 'ssBigBAF'),
#-----------------------------------------

function(object,
         n = c(10),
         mcSamples = 1,                   #number of MC samples
         alpha = 0.05,                    #two-tailed alpha level
         replace = TRUE,                  #MC sample with[out] replacement
         boot = FALSE,
         numBSS = 100,                    #number of bootstrap samples
         description = 'Big BAF Sampling Monte Carlo Simulation',
         startSeed = 123,
         runQuiet = FALSE,
         sshh = FALSE, 
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Comparisons with the variance estimates to the population (sampling surface)
#   variance are not useful. Rather, in the MC setting, we resort to the usual
#   % of times the confidence intervals caught the population mean stats
#   for each of the variances.
#
#***Note that we must catch samples that were generated where only one point
#   or no points fall within inclusion zones. This can easily happen with a
#   sparse population, and/or with small sample sizes. The three cases are...
#
#   (1) zero: generates NaNs in variances from divide by zero
#   (2) one:  generates NAs due to no variance for vbars, delta & Goodman
#   (3) two or more: but all within the same tree so no variation; this will
#       give a zero for variance of vbars and close to that for (8.42)
#
#   Arguments...
#     object = an object of class "ssBigBAF"
#     n = a vector of sample sizes (points/cells)
#     mcSamples = number of Monte Carlo reps
#     alpha = the two-tailed alpha level
#     replace = see ?sample; sample with or w/o replacement from the surface
#     boot = TRUE: include bootstrap & jakknife intervals; FALSE: don't
#     numBSS = number of bootstrap samples per MC+sample size draw
#     description = some text
#     startSeed = see initaRandomSeed()
#     runQuiet = TRUE: lots of feedback; FALSE: none
#     sshh = FALSE: a little feedback even if runQuiet=TRUE, otherwise nothing extra
#     ... = gobbled
#
#   Returns...
#
#     An object of class "monteBigBAF" invisibly.
#
#   The "monte" code parallels that in sampSurf::monteNTSample code; but
#   we are using more than one variable and multiple variances for volume
#   so lists of variances are used here.
#
#   Resampling methods are performed if desired using monteBootBB(), which calls
#   Efron's bcaboot() routine. See the caveat in the following...
#
#***Important: Note that bcajack() currently does not seem to return the jackknife
#              mean, so I am calculating it. I could simply change their routine
#              and this might be more efficeint in the end.
#
#Author...									Date: 14-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   required packages...
#
#    require(sampSurf, quietly = TRUE)
#    if(boot)
#      require(bcaboot, quietly = TRUE)
 
#
#   a few checks...
#
    #checks from monteNTSample...  
    nn = length(n)
    if(nn < 1 || !is.numeric(n))
      stop('Please specify \"n\" as an integer vector!')
    n = sort(round(n))                        #integer sample sizes only, sort drops NAs by default
    n.names = paste('n', n, sep='.')                  #to be compatible with sampSurf monte stuff
    names(n) = n.names
    if(mcSamples < 2 || !is.numeric(mcSamples) || length(mcSamples)>1) #may eventually allow vector
      stop('Please specify \"mcSamples\" as an integer scalar >= 2!')


#
#   initialize the random number generator...
#
    ranSeed = initRandomSeed(startSeed)

#
#   unpack the object for simpler use...
#
    ss.bb.vol = object@ss.bb.vol
    ss.bb.ba = object@ss.bb.ba
    ss.ct.ba = object@ss.ct.ba
    ss.ct.vol = object@ss.ct.vol
    
    estimate = ss.bb.vol@estimate                               #volume in BB

    baf.c = ss.ct.ba@izContainer@iZones[[1]]@angleGauge@baf
    baf.v = ss.bb.ba@izContainer@iZones[[1]]@angleGauge@baf


#
#   get the list of cells within stems for the desired surfaces; remember that
#   background cells are not included in either list...
#
    csl.bb = object@csl.bb                              #BB "ssCellStemList" object
    csl.ct = object@csl.ct                              #Count "ssCellStemList" object

    
#
#----------------------------------------------------------------------------------------------
#           Population Stuff...                                                                
#----------------------------------------------------------------------------------------------
#
#   extract the trees from the sampling surface objects...
#
    gt = getTreesBB(ss.bb.vol)                             #get various tree structures...
    trees = gt$trees                                       #data frame with "vol" & "vbar"
    strees = gt$strees                                     #standingTrees container object
    nTrees = gt$nTrees                                     #number of trees
    tree.ids = trees$id                                    #tree ids


#
#   note that the mean ba from the count surface is over all cells, including
#   background, same as population surface volume...
#
    N = ncell(ss.bb.vol@tract)                           #total number of cells/points
    if(any(n >= N))
      stop('Sample sizes \"n\" can not be larger than the population size \"N\"!')
    ba.count = ss.ct.ba@surfStats$mean                   #basal area from count surface
    var.count = ss.bb.ba@surfStats$var/N                 #count surface ba variance of the mean
    popVolume =  strees@stats['total','volume']          #volume of all trees
    popVolume.var = ss.bb.vol@surfStats$var/N            #volume variance of the mean

#    popVbar.bb = mean(getValues(ss.bb.vol@tract))/mean(getValues(ss.bb.ba@tract)) #from surface
    popVbar.trees = popVolume/strees@stats['total','basalArea']                   #from trees
    
#    pop.tvbs = getValues(ss.bb.vol@tract)/getValues(ss.bb.ba@tract)
#    idx.tvbs = ifelse(is.nan(pop.tvbs), FALSE, TRUE)                              #zero/zero check
#    popCorr.bb = cor(pop.tvbs[idx.tvbs], getValues(ss.bb.ba@tract)[idx.tvbs])     #from surface
    popCorr.trees = with(trees, cor(vbar, ba))                                    #from trees

#
#   ***these are just calculated, I have not used them below yet*********
#   the reson for this is that it is not clear whether they should also
#   be applied to the product variance -- probably not, though they are
#   also variances of the mean. More importantly, the sampling fraction
#   n/N is normally going to be smaller than 5-10%, and as noted in Cochran
#   p. 25, it can be ignored...
#
    fpc = (N - n)/N
    names(fpc) = n.names

#
#   t-values for alpha level desired...
#
    alpha2 = alpha/2
    t.values = qt(1-alpha2, n - 1)
    names(t.values) = n.names

#
#   all four SS will have the same extents, & so share the same cell numbers; we will
#   sample out of the cell numbers and index the tract values from these; the following
#   is a bit of overkill since the cells are stored sequentially from the top left
#   by row to the bottom right in raster, but it is probably best to be careful and
#   rely on raster to do this...
#
    cellNums = cellsFromExtent(ss.bb.vol@tract, extent(ss.bb.vol@tract))



#
#----------------------------------------------------------------------------------------------
#   This is where the MC loop resides...                                                       
#----------------------------------------------------------------------------------------------
#
#   first, set up the storage for the results--mostly lists of data frames...
#
    mc.names = paste('mc', 1:mcSamples, sep = '.')  #names for the rows...
    rvec = data.frame(matrix(NA_real_, nrow = mcSamples, ncol = nn))
    colnames(rvec) = n.names
    rownames(rvec) = mc.names
    #the means useful in estimation...
    means = list(vol.bb = rvec,        #the big baf sample mean for volume
                 ba.bb = rvec,         #the big baf sample mean for basal area
                 ba.ct = rvec,         #the small/count baf sample mean for basal area 
                 vol.ct = rvec,        #the small/count baf sample mean for volume
                 tvbar = rvec          #the mean tree vbar for "in" trees on big BAF points
                )
    #variance and stDevs of the actual sample observations...
    vars = means
    stDevs = vars
    #the estimated variances of the mean corresponding to the big BAF volume estimate...
    varMeans = list(delta = rvec,      #the Delta Method
                    goodm = rvec,      #Goodman's method
                    delta.r = rvec,    #Delta using the ratio variance of vbars     ****
                    goodm.r = rvec,    #Goodman's using the ratio variance of vbars ****
                    vol.bb = rvec,     #H-T volume variance on BB points
                    vol.ct = rvec      #H-T volume variance on count points
                    #Delta = rvec #******with covariance
                   )
    stErrs = varMeans
    #associated lower and upper confidence intervals limits...
    lowerCIs = varMeans
    upperCIs = varMeans
    #whether the CIs caught the population mean (TRUE) or no (FALSE)...
    lvec = rvec
    lvec[] = FALSE
    caught = list(delta = lvec,         #all as in varMeans, but with logical seed
                  goodm = lvec, 
                  delta.r = lvec,    #Delta using the ratio variance of vbars     ****
                  goodm.r = lvec,    #Goodman's using the ratio variance of vbars ****
                  vol.bb = lvec, 
                  vol.ct = lvec
                 )
    ##******caught = list(delta = lvec, goodm = lvec, vol.bb = lvec, vol.ct = lvec, Delta = lvec) #******
    #other variances of the mean of possible interest...
    otherVarms = list(ba.bb = rvec,         #variance of BB basal area
                      ba.ct = rvec,         #variance of count ba
                      tvbar = rvec,         #variance of tree vbars
                      tvbar.r = rvec,       #ratio variance of vbars      ****
                      dm.tvbar = rvec,      #delta method tree vbar volume component
                      dm.tvbar.r = rvec,    #delta method tree vbar volume component
                      dm.ba = rvec          #delta method point basal area volume component
                     )    
    #the number of vbar sample trees selected in a given sample...
    n.tvbar = list(count = rvec,   #for the smaller/count BAF
                   bb = rvec       #for the big BAF
                  )
    #the covariance and correlation between the vbar & ba on the BB points...
    covs = list()            #placeholder for PBDM covariances (added 5-May-2020)
    corrs = list(tvbar.ba = rvec,          #individual tree vbars and tree basal areamc.samples
                 Tvbar.ba = rvec,          #aggregate tree vbar & ba for each point
                 pvbar.ba = rvec,          #cell/point aggregate vbars & total cell ba
                 mvbar.ba = rvec,
                 Pvbar.ba = rvec,
                 Mvbar.ba = rvec
                )
#
#   augment some of the above for bootstrapping and jackknifing...
#
    if(boot) {
      means = c(means, list(boot = rvec, jack = rvec))
      varMeans = c(varMeans, list(boot = rvec, jack = rvec))
      stErrs = c(stErrs, list(boot = rvec, jack = rvec))
      lowerCIs = c(lowerCIs, list(boot = rvec, jack = rvec))
      upperCIs = c(upperCIs, list(boot = rvec, jack = rvec))
      caught = c(caught, list(boot = lvec, jack = lvec))
    }


#
#   this will hold the actual cell numbers used in each individual MC sample by
#   sample size...
#
    mc.samples = vector('list', nn)
    names(mc.samples) = n.names
    for(i in seq_len(nn)) {
      df = data.frame(matrix(NA_real_, nrow = n[i], ncol = mcSamples))
      #rownames(df) = paste('c', 1:mcSamples, sep = '.')
      colnames(df) = mc.names
      mc.samples[[i]] = df
    }
    rm(rvec, lvec, df)


#-------------------------------------------------------------------------------
#   Monte Carlo and sample size loops; note that when the MC loop is the inner
#   different runs with, say, one less sample size are comparable...
#
    if(!sshh) cat('\nn=')
    for(j in seq_len(nn)) {                               #sample size first
      n.j = n[j]
      if(!sshh)
        cat(n.j,', ',sep='')
      for(i in seq_len(mcSamples)) {                      #MC reps w/in sample size
 
        #sample the cell numbers to use as indexes; these are the point samples...
        sdx = sample(cellNums, n.j, replace = replace)    #obviously, this can contain background cells

        #summarize for each surface...
        stats = monteStatsBB(sdx, object)
        stats.ct = stats$stats.ct                         #for count baf
        stats.bb = stats$stats.bb                         #for big baf
        mc.samples[[j]][,i] = stats$samples$sdx           #sorted cell numbers
        

        #get the weighted tree vbar stats...
        mts.bb = monteTreeStatsBB(sdx, object, 'bigBAF')                     #for BB cells
        #tree.vbars = with(mts.bb$df.trees, vbar[idx]) #not used!
        tvbar.mean = mts.bb$treeStats['mean', 'vbar']
        tvbar.var = mts.bb$treeStats['var', 'vbar']
        tvbar.varm = mts.bb$treeStats['varm', 'vbar']
        vars$tvbar[i,j] = tvbar.var
        stDevs$tvbar[i,j] = sqrt(tvbar.var)
        #number of trees selected including multiple sample trees on different points...
        mts.ct = monteTreeStatsBB(sdx, object, 'count')                     #Count sample
        n.tvbar$count[i,j] = mts.ct$n.v
        n.tvbar$bb[i,j] = mts.bb$n.v                                      #BB sample

        #means...
        means$vol.bb[i,j] = stats.ct['mean', 'ba']*tvbar.mean    #bb product mean, not surface vol
        means$ba.bb[i,j] = stats.bb['mean', 'ba']
        means$vol.ct[i,j] = stats.ct['mean', 'vol']
        means$ba.ct[i,j] = stats.ct['mean', 'ba']
        means$tvbar[i,j] = tvbar.mean
        
        #observed sample variances and standard deviations...
        vars$vol.bb[i,j] = stats.bb['var', 'vol']                 #of surface vol, not from bb
        vars$ba.bb[i,j] = stats.bb['var', 'ba']
        vars$vol.ct[i,j] = stats.ct['var', 'vol']
        vars$ba.ct[i,j] = stats.ct['var', 'ba']
        stDevs$vol.bb[i,j] = sqrt(vars$vol.bb[i,j])
        stDevs$ba.bb[i,j] = sqrt(vars$ba.bb[i,j])
        stDevs$vol.ct[i,j] = sqrt(vars$vol.ct[i,j])
        stDevs$ba.ct[i,j] = sqrt(vars$ba.ct[i,j])

        #variances of the mean & standard errors on volume from bb...
        #(note the NA catch for 0 or 1 vbar trees below)
        dm.tvbar = stats.ct['mean', 'ba']^2*tvbar.varm         #volume component from vbars
        dm.ba = tvbar.mean^2*stats.ct['varm', 'ba']            #volume component from basal areas
        if(!is.na(dm.tvbar))                                   #will be NA if 0 or 1 vbar trees  *++*
          vt.delta = dm.tvbar + dm.ba
        else                                                   #this will keep it all from being NA  *++*
          vt.delta = dm.ba                                     #i.e., just use this component  *++*
        varMeans$delta[i,j] = vt.delta
        if(!is.na(dm.tvbar))                                   #will be NA if 0 or 1 vbar trees  *++*
          varMeans$goodm[i,j] = vt.delta - stats.ct['varm', 'ba']*tvbar.varm
        else                                                   #this will keep it all from being NA  *++*
          varMeans$goodm[i,j] = vt.delta                       #i.e., just use this component  *++*
        stErrs$delta[i,j] = sqrt(vt.delta)
        stErrs$goodm[i,j] = sqrt(varMeans$goodm[i,j])
        #as above, but using the ratio variance estimator for tree vbars...  ****
        #(ratio is okay if one vbar tree, and will be 0.0 for 0 vbar trees unlike the above)
        dm.tvbar.r = stats.ct['mean', 'ba']^2*stats.bb['varm', 'vbar']
        vt.delta.r = dm.tvbar.r + dm.ba                        #of course, ba remains as above
        varMeans$delta.r[i,j] = vt.delta.r
        varMeans$goodm.r[i,j] = vt.delta.r - stats.ct['varm', 'ba']*stats.bb['varm', 'vbar']
        stErrs$delta.r[i,j] = sqrt(vt.delta.r)
        stErrs$goodm.r[i,j] = sqrt(varMeans$goodm.r[i,j])
        
        #H-T variance G&V (8.31) for bb...
        varMeans$vol.bb[i,j] = stats.bb['varm', 'vol']
        stErrs$vol.bb[i,j] = sqrt(varMeans$vol.bb[i,j])
        #H-T variance G&V (8.31) for count...
        varMeans$vol.ct[i,j] = stats.ct['varm', 'vol']
        stErrs$vol.ct[i,j] = sqrt(varMeans$vol.ct[i,j])
        
        #bootstrapping and jackknifing--means and variances...
        if(boot) {
          capture.output({ #bcajack() prints some garbage that we don't need to see...
            mcb = monteBootBB(stats$samples, numBSS, alpha)
          })
          means$boot[i,j] = mcb$stats['est','theta']
          varMeans$boot[i,j] = mcb$stats['est','sdboot']*mcb$stats['est','sdboot']
          stErrs$boot[i,j] = mcb$stats['est','sdboot']
          means$jack[i,j] = mcb$jackmean  
          varMeans$jack[i,j] = mcb$stats['est','sdjack']*mcb$stats['est','sdjack']
          stErrs$jack[i,j] = mcb$stats['est','sdjack']
        }
        
        #------------------------------------------------------
        #CIs and catch for each variance method on bb volume...
        #------------------------------------------------------
        #for the delta method...
        x = .ciCaught(means$vol.bb[i,j], stErrs$delta[i,j], t.values[j], popVolume)
        lowerCIs$delta[i,j] = x$lci
        upperCIs$delta[i,j] = x$uci
        caught$delta[i,j] = x$catch
        #for Goodman's method...
        x = .ciCaught(means$vol.bb[i,j], stErrs$goodm[i,j], t.values[j], popVolume)
        lowerCIs$goodm[i,j] = x$lci
        upperCIs$goodm[i,j] = x$uci
        caught$goodm[i,j] = x$catch
        #for the delta method with ratio variance... ****
        x = .ciCaught(means$vol.bb[i,j], stErrs$delta.r[i,j], t.values[j], popVolume)
        lowerCIs$delta.r[i,j] = x$lci
        upperCIs$delta.r[i,j] = x$uci
        caught$delta.r[i,j] = x$catch
        #for Goodman's method with ratio variance...  ****
        x = .ciCaught(means$vol.bb[i,j], stErrs$goodm.r[i,j], t.values[j], popVolume)
        lowerCIs$goodm.r[i,j] = x$lci
        upperCIs$goodm.r[i,j] = x$uci
        caught$goodm.r[i,j] = x$catch
        #for H-T method on the large BAF...
        x = .ciCaught(means$vol.bb[i,j], stErrs$vol.bb[i,j], t.values[j], popVolume)
        lowerCIs$vol.bb[i,j] = x$lci
        upperCIs$vol.bb[i,j] = x$uci
        caught$vol.bb[i,j] = x$catch
        #for H-T method on the small/count BAF...
        x = .ciCaught(means$vol.ct[i,j], stErrs$vol.ct[i,j], t.values[j], popVolume)
        lowerCIs$vol.ct[i,j] = x$lci
        upperCIs$vol.ct[i,j] = x$uci
        caught$vol.ct[i,j] = x$catch
        
        #bootstrapping and jackknifing...
        if(boot) {
          alpha.lo = as.character(alpha2)
          alpha.hi = as.character(1 - alpha2)
          lowerCIs$boot[i,j] = mcb$lims[alpha.lo, 'bca']
          upperCIs$boot[i,j] = mcb$lims[alpha.hi, 'bca']
          if(!is.na(lowerCIs$boot[i,j]) && !is.na(upperCIs$boot[i,j])) {
            if(lowerCIs$boot[i,j] <= popVolume && popVolume <= upperCIs$boot[i,j])
              caught$boot[i,j] = TRUE
          } #NA if
          #determine standard error-based jackknife intervals too...
          x = .ciCaught(means$jack[i,j], stErrs$jack[i,j], t.values[j], popVolume)
          lowerCIs$jack[i,j] = x$lci
          upperCIs$jack[i,j] = x$uci
          caught$jack[i,j] = x$catch
        } #boot
        
        #other variances of the mean...
        otherVarms$ba.bb[i,j] = stats.bb['varm', 'ba']
        otherVarms$ba.ct[i,j] = stats.ct['varm', 'ba']
        otherVarms$tvbar[i,j] = tvbar.varm
        otherVarms$tvbar.r[i,j] = stats.bb['varm', 'vbar']    #ratio estimate of tree vbar varm ****
        #components from the delta method...
        otherVarms$dm.tvbar[i,j] = dm.tvbar 
        otherVarms$dm.tvbar.r[i,j] = dm.tvbar.r
        otherVarms$dm.ba[i,j] = dm.ba
        

        #------------------------
        #correlation section...  
        #------------------------
        #total individual tree vbars by cell/point & expanded total ba per cell/point;
        #these are both totals, i.e., sums for a given point so they are the point-wise
        #averages for a single point...
        corrs$pvbar.ba[i,j] = cor(mts.bb$df.cells$vbar, stats$samples$ba.cts)
        corrs$Pvbar.ba[i,j] = cor(mts.bb$df.cells$vbar, stats$samples$ba.bbs)    #BB over all points
        ##corrs$pvbar.ba[i,j] = cor(mts.bb$df.cells$vbar, mts.ct$df.cells$ba) #no: both tree-based
        
        #individual tree vbars & ba, over each tree in the sample, not by cell/point...
        #because there are different no of trees sampled per point with the BB and CT
        #samples, we have to correlate over just, e.g., BB...
        corrs$tvbar.ba[i,j] = mts.bb$treeStats['corr', 'vbar']
        #now aggregate tree vbar and ba point-wise rather than tree-wise...
        corrs$Tvbar.ba[i,j] = mts.bb$cellStats['corr', 'vbar']
        
        #the following is like pvbar, but for the mean tree per point: the mean tree
        #vbar and the mean tree total basal area (i.e., it reduces to the baf regardless
        #of the number of tree sampled per point)...
        m.vb.bb = with(mts.bb$df.cells, ifelse(is.nan(vbar/no.trees), 0.0, vbar/no.trees))
        #the following with() is to get the no.trees required...
        m.ba.ct = with(mts.ct$df.cells, ifelse(is.nan(stats$samples$ba.cts/no.trees), 0.0, 
                                                      stats$samples$ba.cts/no.trees)) #constant baf
        corrs$mvbar.ba[i,j] = cor(m.vb.bb, m.ba.ct)
        #BB means over all points, including zero-count points...
        M.vb.bb = with(mts.bb$df.cells, ifelse(is.nan(vbar/no.trees), 0.0, vbar/no.trees))
        #the following with() is to get the no.trees required...
        M.ba.bb = with(mts.bb$df.cells, ifelse(is.nan(stats$samples$ba.bbs/no.trees), 0.0, 
                                                      stats$samples$ba.bbs/no.trees)) #constant baf
        corrs$Mvbar.ba[i,j] = cor(M.vb.bb, M.ba.bb)
        
      } #sample size loop
    } #MC loop  


#
#   average summary statistics--grand means over each sample for each n;
#   note that it also can contain jackknife & bootstrap results...
#
    gm.means = .grandMeans(means, n.names)
    gm.vars = .grandMeans(vars, n.names)
    gm.stDevs = .grandMeans(stDevs, n.names)
    gm.varMeans = .grandMeans(varMeans, n.names)
    gm.stErrs = .grandMeans(stErrs, n.names)
    gm.lowerCIs = .grandMeans(lowerCIs, n.names)
    gm.upperCIs = .grandMeans(upperCIs, n.names)
    gm.caught = .grandMeans(caught, n.names) * 100
    gm.otherVarms = .grandMeans(otherVarms, n.names)
    gm.corrs = .grandMeans(corrs, n.names)
    gm.all = list(gm.means = gm.means, gm.vars = gm.vars, gm.stDevs = gm.stDevs,
                  gm.varMeans = gm.varMeans, gm.stErrs = gm.stErrs,
                  gm.lowerCIs = gm.lowerCIs, gm.upperCIs = gm.upperCIs,
                  gm.caught = gm.caught, gm.otherVarms = gm.otherVarms,
                  gm.corrs = gm.corrs
                 )

#
#   create the object (not quite complete at this point: will fail validObject)...
#
    ssBB.mc = new('monteBigBAF',
                  description = description,
                  estimate = estimate,
                  mcSamples = mcSamples,
                  n = n,
                  fpc = fpc,
                  alpha = alpha,
                  ranSeed = ranSeed,
                  t.values = t.values,
                  boot = boot,
                  numBSS = numBSS,
                  replace = replace,
                  means = means,
                  vars = vars,
                  stDevs = stDevs,
                  varMeans = varMeans,
                  stErrs = stErrs,
                  lowerCIs = lowerCIs,
                  upperCIs = upperCIs,
                  caught = caught,
                  n.tvbar = n.tvbar,
                  otherVarms = otherVarms,
                  corrs = corrs,
                  gm.all = gm.all,
                  mc.samples = mc.samples
                 ) #new

#
#   the approximate sampling variances of the means from the MC replications;
#   computed for each of the four sampSurf objects...
#
    sm.all = .samplingVarMeans(ssBB.mc) 
    ssBB.mc@sm.all = sm.all   

#
#   augment the object with the PBDM statistics...
#
    ssBB.mc = monte(ssBB.mc, object, runQuiet = TRUE)

#
#   finally print the population totals if desired...
#
    if(!runQuiet) {
      cat('\nPopulation of', nTrees, 'trees...')
      cat('\nPopulation size N =', N, '(number of cells)')
      cat('\nSample size =', n)
      cat('\nSampling fraction n/N =', n/N)
      cat('\nSamples drawn', ifelse(replace, 'with', 'without'), 'replacement')
      cat('\nBasal Area...')
      cat("\n  Population tree ba =", ss.bb.ba@surfStats$popTotal)  #==sum(trees$ba)
      cat("\n  mean ss.bb.ba =", ss.bb.ba@surfStats$mean)
      cat("\n  mean ss.ct.ba =", ba.count)
      cat('\nVolume...')
      cat('\n  Population tree volume =', popVolume)
      cat("\n  mean ss.bb.vol volume =", ss.bb.vol@surfStats$mean)
      cat('\n    %diff pop vol =',  pctDiff(ss.bb.vol@surfStats$mean, popVolume))
      cat('\n  Population mean Vbar =', popVbar.trees)
      cat('\nCorrelation (tree vbar & ba)...')
      cat('\n  Population on individual trees =', popCorr.trees)
#      cat('\n  BB surface =', popCorr.bb)
      cat('\n')
    }
   if(!sshh) 
     cat('\n')
             
    return(invisible(ssBB.mc))

}   #monte constructor for 'ssBigBAF' objects
)   #setMethod
    

