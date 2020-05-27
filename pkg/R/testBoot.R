testBoot = function(n = 1000,
                    B = 1000,
                    alpha = 0.05,
                    mbm = FALSE,
                    parallel = 'no',
                    ncpus = 3L,
                    times = 10L,
                    unit = 's',
                    runQuiet = FALSE,
                    startSeed = 245,
                    ...
                   )
{
#---------------------------------------------------------------------------
#
#   This will just check the two packages (bcaboot, and boot) for the
#   calculation of "bca" confidence intervals. Interest lies in both
#   how close the estimated intevals are in both cases, and also the time
#   it takes to calculate the intervals. Note below that for comparisons
#   the random number generator is reset using the same seed before each
#   of the bootstrap iterations.
#
#   Note that I assume bcajack is at a slight disadvantage from the start
#   since it is also calculating the jackknife standard error. But, boot
#   takes two function calls, one to bootstrap and one to calculate the
#   intervals. So perhaps they both will be close.
#
#   Arguments...
#     n = size for the sample that we will bootstrap
#     B = number of bootstrap sample replicates
#     alpha = two-tail level
#     mbm = TRUE: run microbenchmark for times; FALSE; no timing
#     parallel = see boot, it does not seem to work at all
#     ncpus = see boot
#     times = the number of microbenchmark replicates
#     units = the units for the report from microbenchmark; 's' = seconds
#     runQuiet = TRUE: sssshhh; FALSE: print results
#     startSeed = the starting seed
#     ... = other arguments passed to boot()
#
#   Returns...
#     a list invisibly with...
#       -- the sample drawn
#       -- results from bcajack()
#       -- results from boot()
#       -- results from boot.ci()
#       -- microbenchmark results for bcaboot()
#       -- microbenchmark results for boot() & boot.ci()
#       -- data frame with the normal theory results from the sample
#       -- data frame with the bootstrap results
#       -- call is from match.call()
#
#   If we set "mbm=TRUE" then a time comparison is made with microbenchmark.
#
#   It is disappointing how different the two sets of bca intervals seem
#   to be. They do "converge" as B gets large, but it must be in the ten
#   to 20 thousand range. The biggest offender is the upper CI difference.
#
#   Here are some results using the default startSeed above...
#
#                              bcajack                  boot       
#       n        B       lower      upper          lower     upper 
#     1000     1000  -0.10200987  0.015298104  -0.098467178 0.027124759
#     1000     2000  -0.10273359  0.017794756  -0.10158977  0.023272337
#     1000    10000  -0.10273359  0.019465835  -0.10256112  0.020762515
#     1000    20000  -0.10284077  0.020526647  -0.10198038  0.02077783 
#     1000    40000  -0.10206452  0.020885554  -0.10038422  0.021616858
#   ...
#     normal theory: -0.10192412  0.02085006
#
#   The means do agree precisely between the two.
#
#   Also, the timing stats seem to prefer boot over bcajack. But another
#   disappointment is that the "parallel" argument in boot does not seem
#   to matter--one gets the same results with or without it set to
#   "multicore" on average.
#
#   Be careful with boot.ci as it requires B>n evidently or you will get
#   an error. One time I got an error because B was not a multiple of n.
#   It seems that making B = 2n is the usual trick; e.g.,...
#
#   https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html
#
#
#Author...									Date: 28-Mar-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#    require(boot, quietly = TRUE)
#    require(bcaboot, quietly = TRUE)
#    if(mbm)
#      require(microbenchmark, quietly = TRUE)
#    if(parallel != 'no')
#      require(parallel, quietly = TRUE)
    
#
#   initialize the random number generator; note that this is also reset for
#   each bootstrap simulation in the sections below...
#
    ranSeed = initRandomSeed(startSeed)     #or simply use set.seed()

#
#   now the sample from which we will bootstrap; calculate the normal theory intervals
#   for comparison...
#
    samp = rnorm(n)
    smean = mean(samp)
    stErr = sqrt(var(samp)/n)
    t.val = qt(1-alpha/2, n-1)
    CIs = c(smean - t.val*stErr, smean + t.val*stErr)
    if(!runQuiet) {
      cat('\nOriginal sample normal theory intervals...')
      cat('\n  Mean =', smean, 'for n =', n)
      cat('\n ', 1-alpha, 'CIs =', CIs)
    }
    df.samp = data.frame(smean, stErr, CIs[1], CIs[2])
    colnames(df.samp) = c('mean', 'stErr', 'ci.lo', 'ci.hi')


    df = data.frame(matrix(NA_real_, nrow = 1, ncol = 6))
    colnames(df) = c('n', 'B', 'bcajack.lo', 'bcajack.hi', 'boot.lo', 'boot.hi')
    df[1,'n'] = n
    df[1,'B'] = B

#
#-------------------------------------------------
#   bcaboot test; en==Efron & Narasimhan...
#-------------------------------------------------
#
    ranSeed2 = initRandomSeed(startSeed) #re-initialize
    alphas = c(alpha/2, 1-alpha/2)
    capture.output({ #bcajack() prints some garbage that we don't need to see...
       res = bcajack(samp, B, mean, alpha = alphas)
    })
    if(!runQuiet) {
      cat('\nbcajack...')
      cat('\n  mean =', res$stats['est',1])
      cat('\n ', 1-alpha, 'bca cis = ', res$lims[as.character(alphas), 'bca'] )
    }
    df[1, c('bcajack.lo', 'bcajack.hi')] = res$lims[as.character(alphas), 'bca']
    
#
#   benchmark bcajack...
#
    if(mbm) {
      capture.output({ 
        mbm.en = microbenchmark({
          res2 = bcajack(samp, B, mean, alpha = alphas)
        }, times = times, unit = unit )
      })
      s = summary(mbm.en, unit=unit)
      if(!runQuiet)
        cat("\n  Unit: ", attr(s, "unit"), "\n", sep = "")
      print(s[,-1])
    }
    else
      mbm.en = NULL


#
#-------------------------------------------------
#   boot; cr==Canty & Ripley...
#-------------------------------------------------
#
    ranSeed3 = initRandomSeed(startSeed) #re-initialize again
    meanboot = function(x,idx) mean(x[idx])  #calculate the mean for each BS sample
    boot.samp = boot(samp, meanboot, B, parallel = parallel, ...)
    boot.mean = mean(boot.samp$t0)
    boot.cis = boot.ci(boot.samp, 1-alpha, type='bca')
    lower = boot.cis$bca[1, 4]
    upper = boot.cis$bca[1, 5]
    if(!runQuiet) {
      cat('\nboot...')
      cat('\n  mean =', boot.mean)
      cat('\n ', 1-alpha, 'bca cis = ', c(lower, upper))
    }
    df[1, c('boot.lo', 'boot.hi')] = c(lower, upper)
    
#
#   benchmark boot...
#
    if(mbm) {
      capture.output({ 
        mbm.cr = microbenchmark({
           boot.samp = boot(samp, meanboot, B, parallel = parallel, ncpus = ncpus, ...)
           boot.cis = boot.ci(boot.samp, 1-alpha, type='bca')
        }, times = times, unit = unit )
               
      })
      s = summary(mbm.cr, unit=unit)
      if(!runQuiet)
        cat("\n  Unit: ", attr(s, "unit"), "\n", sep = "")
      print(s[,-1])
    }
    else
      mbm.cr = NULL


       
    if(!runQuiet)
      cat('\n')  
    return(invisible(list(samp = samp,
                          res = res,
                          boot.samp = boot.samp,
                          boot.cis = boot.cis,
                          mbm.en = mbm.en,
                          mbm.cr = mbm.cr,
                          df.samp = df.samp,
                          df = df,
                          call = match.call()
                         )
                     )
          )
    
} #testBoot

