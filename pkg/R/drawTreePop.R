drawTreePop = function(tract,        #units of tract determines the output units
                       solidTypes = c(1.5, 3),
                       topDiams = c(0,0),
                       B = 80,       #80 ft^2 provides a reasonable number of trees
                       hgt.sd = 6,   #feet
                       #Weibull parameters...
                       a = 4.0, #inches
                       b = 8.0, #inches
                       c = 2.0,
                       #spatial package...
                       inhibitDist = 3,     #in unitsOut == tract
                       #other
                       showPlot = TRUE,
                       startSeed = 144,
                       runQuiet = FALSE,
                       ...
                      )
{
#---------------------------------------------------------------------------
#
#   This routine will set up a reasonable simulated population of trees
#   for northern hardwoods based on the Weibull diameter distribution 
#   parameters and height equation from Fast and Ducey (2011). It is based on
#   the shortleaf pine version that was used for the Lynch & Gove Antithetic
#   CHS paper, CJFR (2013).
#
#   The spatial coordinates are drawn within the tract using an inhibition
#   process; the "spatial" package code is used here because it is faster than
#   the corresponding routine in the spatstats package.
#
#   The input units are always "English" for everything but the tract, which can
#   be either; the output units are the same as that of the tract. Therefore,
#   if the tract is metric, then the tree list is also. 
#
#   It is used to set up the synthetic population ...
#
#   Arguments...
#     tract = a "Tract" object in "English" or 'metric'
#     solidTypes, topDiams = see sampleTree() help
#     B = in English, B is BPA 
#     hgt.sd = standard deviation for height for rnorm(0,hgt.sd) perturbations
#     a,b,c = Weibull location, scale and shape parameters ("English")
#     inhibitDist = see ?SSI, inhibition distance (min dist between trees)
#     showPlot = TRUE: plot the sampled dbh distn against the Weibull;
#                FALSE: no plot
#     startSeed = see ?initRandomSeed
#     runQuiet = TRUE: shhhh; FALSE: !shhhh
#     ... = passed to sampleTrees() and hist()
#
#   Returns (invisibly)...
#     A data frame with the synthetic trees that can be cast to standingTrees
#
#Author...									Date: 28-Mar-2017
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
    #require(sampSurf, quietly = TRUE)
#    require(spatial, quietly = TRUE)

#
#   some useful conversions...
#
    acre2ha = 2.4710439
    in2cm = 2.54
    ft2m = 0.3048
    sfa2smh = 0.22956841   #1/4.356
    smpha = 10000
    sfpac = 43560

    unitsOut = tract@units
    if(unitsOut == 'metric')
      hgt.units = 'meters'
    else
      hgt.units = 'feet'

#
#   the tract will always be in the desired output units; if it is buffered,
#   used that as the bbox; note that for a buffered tract, the tractArea is
#   that area INSIDE the buffer, since it is used to draw the correct size
#   for the tree population...
#
    if(is(tract, 'bufferedTract')) {
      buff = TRUE
      bbox = tract@bufferRect
      bboxArea = (bbox[1,2] - bbox[1,1])^2
      tractArea = bboxArea / ifelse(unitsOut == 'metric', smpha, sfpac)
    }
    else {
      buff = FALSE
      bbox = bbox(tract)
      tractArea = area(tract) / ifelse(unitsOut == 'metric', smpha, sfpac)
    }
      

    initRandomSeed(startSeed)

#
#   all equations were fitted in English units...
#
    K = pi/(4*144)
    mu2 = b^2*gamma(2/c+1) + 2*b*a*gamma(1/c+1) + a^2   #qmsd^2
    N = round(B/(K*mu2))
    qmsd = sqrt(mu2)
    if(!runQuiet) {
      cat('\nInput specs...')
      cat('\n  Number of trees/acre =',N)
      cat('\n  Quadratic msd =', qmsd, 'in')
      cat('\n  Basal area/acre =', B)
      cat('\n  Weibull shape =', c)
      cat('\n  Weibull scale =', b, 'in')
      cat('\n  Weibull shift =', a, 'in')
      cat('\nOutput specs...')
      cat('\n  units =', unitsOut)
      if(unitsOut == 'metric') {
        N = round(N*acre2ha)
        B = B*sfa2smh
        qmsd = qmsd*in2cm
        cat('\n  Number of trees/ha =',N)
        cat('\n  Quadratic msd =', qmsd, 'cm')
        cat('\n  Basal area/ha =', B)
        cat('\n  Weibull scale =', b*in2cm, 'cm')
        cat('\n  Weibull shift =', a*in2cm, 'cm')
      }
    }
    else if(unitsOut == 'metric') {
      N = round(N*acre2ha)
      B = B*sfa2smh
      qmsd = qmsd*in2cm
    }
    
#
#   Now, put required stuff on full tract basis in case it's not one acre (ha)...
#
    N = N * tractArea

    if(!runQuiet) {
      cat('\n\n  Totals for the tract...')
      if(buff)
        cat('\n  --Tract area (inside the buffer) =', tractArea)
      else
        cat('\n  -- Tract area =', tractArea)
      cat('\n  --Total N for above area =', N)
    }

#
#   set up the data frame with random stuff, then replace the columns as needed below...
#
    trees = sampleTrees(N, sampleRect=bbox, solidTypes=solidTypes, topDiams=topDiams, ...)

    dbh = rweibull(N, shape=c, scale=b) + a

    #totHgt = 4.5 + 5.7331*dbh*exp(-0.0217*dbh)       #English coefficients for hemlock
    #trees$species = 'hemlock'
    totHgt = 4.5 + 150.9770/(1 + 1/(0.0827*dbh^0.8660))  #English coefficients for all spp
    if(hgt.sd > 0.0) {
      e = rnorm(N, 0.0, hgt.sd)
      totHgt = totHgt + e
      if(!runQuiet) {
        if(unitsOut == 'metric')
          hgt.mult = ft2m
        else
          hgt.mult = 1.0
        hgtOut = hgt.units
        cat('\n  --Height perturbations with sd =', hgt.mult*hgt.sd, hgt.units, 'added.')
      }
    }
    trees$species = 'NHdwds'

#
#   output units to metric?...
#
    if(unitsOut == 'metric') {
      dbh = dbh*in2cm
      totHgt = totHgt*ft2m
    }
    trees$dbh = dbh
    trees$height = totHgt
    trees$units = unitsOut

#
#   generate spatial locations within the tract bbox...
#
    spatial::ppregion(bbox[1,1],bbox[1,2],bbox[2,1],bbox[2,2]) #set the domain for SSI()
    xy = spatial::SSI(N, inhibitDist)
    trees$x = xy$x
    trees$y = xy$y
    
    
#
#   show a plot with sampled trees vs. Weibull if desired...
#
    if(showPlot) {
      conv = ifelse(unitsOut == 'metric', in2cm, 1.0)
      a.plt = a * conv
      b.plt = b * conv
      d.min = 0.9 * a.plt
      d.max = 1.1 * ceiling(max(dbh))
      hist(dbh, breaks = seq(d.min, d.max, by = 1), col = 'gray90', xlim = c(d.min, d.max),
           ...) 
      x = seq(d.min, d.max, by = 0.5)
      lines(x, dweibull(x - a.plt, c, b.plt)*N, col='red')
    }

    if(!runQuiet) {
      if(unitsOut == 'metric')
        K = pi/(4*10000)
      else
        K = pi/(4*144)
      cat('\n  --Total Basal area sampled =', sum(dbh*dbh*K))
      cat('\n')
    }

    
    return(invisible(trees))
}   #drawTreePop

