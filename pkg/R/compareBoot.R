compareBoot = function(n = c(100, 200),
                       B = c(1000, 1200),
                       tablePath = '',
                       tableName = '',
                       digits = c(0,0,0,rep(3,6)),
                       caption = 'Bootstrap Confidence Interval Comparison',
                       label = 'tab:bootTable',
                       ssshhh = FALSE,
                       ...
                      )
{
#---------------------------------------------------------------------------
#
#   This will run testBoot() on each n+B combination. Thus, it can be used
#   to look at the differences between the two bca routines, bcaboot::bcajack
#   and boot:boot.ci for confidence interval calculation.
#
#   Arguments...
#     n = size for the sample that we will bootstrap
#     B = number of bootstrap sample replicates
#     tablePath = the relative path to tableName
#     tableName = '' for no hardcopy table output; otherwise a *.tex name
#     digits = the number of display digits for the xtable
#     caption = a caption for the xtable
#     label = a LaTeX table label for the xtable
#     ssshhh = runQuiet is set to TRUE in the call to testBoot, so this 
#              is a local version only
#     ... = passed to testBoot()
#
#   Returns...
#     a list invisibly with...
#       -- a data frame with CI results
#       -- and xtable objects with the above in tex format
#
#Author...									Date: 19-Apr-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   set up the data frame for the results...
#
    n.n = length(n)
    n.B = length(B)
    df.names = c('n', 'B', 'Norm.lo', 'Norm.hi', 'bcajack.lo', 'bcajack.hi',
                 'boot.lo', 'boot.hi')
    n.cols = length(df.names)
    df = data.frame(matrix(NA_real_, nrow = n.n*n.B, ncol = n.cols))
    names(df) = df.names    

#
#   loop through all sample sizes and bootstrap sample sizes...
#
    ii = 1                                #index counter
    for(i in seq_len(n.n)) {
      for(j in seq_len(n.B)) {
        z = testBoot(n[i], B[j], runQuiet = TRUE, ...)
        df[ii, 1] = n[i]
        df[ii, 2] = B[j]
        df[ii, 3] = z$df.samp[1, 'ci.lo']
        df[ii, 4] = z$df.samp[1, 'ci.hi']
        df[ii, 5] = z$df[1, 'bcajack.lo']
        df[ii, 6] = z$df[1, 'bcajack.hi']
        df[ii, 7] = z$df[1, 'boot.lo']
        df[ii, 8] = z$df[1, 'boot.hi']
        ii = ii + 1
      } #j
    } #i


#
#   now save it to a tex file if desired; make an xtab data frame if xtable is available...
#
    if(requireNamespace("xtable", quietly=TRUE)) {
      latexColumns = colnames(df)
      xdf = df[, latexColumns]               #desired subset of all columns
      xtab = xtable::xtable(xdf, caption = caption, label = label, 
                            align = rep('c', n.cols+1), digits = digits)
      if(nchar(tableName) > 0) {
        fn = file.path(getwd(), tablePath, tableName)
        if(!ssshhh)
          cat('\nExporting', fn, '\n')
        print(xtab, 'latex', fn, include.rownames = FALSE)
      }
    }
    else
      xtab = NA

    return(invisible(list(df = df,
                          xtab = xtab
                         )
                    )
          )
}   #compareBoot

