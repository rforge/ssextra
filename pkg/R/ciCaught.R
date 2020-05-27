.ciCaught = function(mean, 
                     stErr, 
                     t.val,
                     popVolume
                    )
{
#---------------------------------------------------------------------------
#
#   This little function computes the CIs and capture outcome for a given 
#   Monte Carlo sample, whose mean and standard error of the mean are passed.
#
#   It is very similar to the in-line function in ssExtra::monteBBConstructor
#   except that here we pass the standard error of the mean instead of the
#   variance of the mean.
#
#   This is used in montePBDM, which is getting large, so I have kept it
#   separate to reduce the clutter in the code.
#
#   Arguments...
#     mean = the mean of a single MC sample
#     stErr = the standard error of the mean above
#     t.val = the two-tailed t-value
#     popVolume = the population volume for the catch
#
#   Returns...
#     A list invisibly with...
#       -- lower CI endpoint
#       -- upper CI endpoint
#       -- TRUE: caught the population mean; FALSE: missed it
#
#Author...									Date: 25-Nov-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   very simple...
#
    ts = t.val*stErr
    lci = mean - ts
    uci = mean + ts
    
#
#   see if we caught the population mean...
#
    catch = FALSE
    if(!is.na(lci) && !is.na(uci)) {
      if(lci <= popVolume && popVolume <= uci)
        catch = TRUE
    }
    
    return(invisible(list(lci = lci,
                          uci = uci,
                          catch = catch
                         )
                    )
          )
}   #.ciCaught


