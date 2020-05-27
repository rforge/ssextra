#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   ssCellStemList class...
#
#   The methods include signatures for...
#     Signature: "object"
#     1. "sampSurf":    Takes one sampSurf object
#
#
#Author...									Date: 1-May-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#   generic definition...
#
setGeneric('ssCellStemList',  
           function(ss, ...) standardGeneric('ssCellStemList'),
             signature = c('ss')
            )


          
#================================================================================
#
# 1. Takes one "sampSurf" object...
#
setMethod('ssCellStemList',
          signature(ss = 'sampSurf'), 
function(ss,
         description = 'list of cells in stems & stems in cells',
         ...
        )
{
#---------------------------------------------------------------------------
#
#   This is the main way to create a simple "ssCellStemList" object...
#
#Author...									Date: 1-May-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   use the basic function to develop the list structures...
#
    lcs = listCellStem(ss, cellsOnly = FALSE)

#
#   get the data frame of the number of cells containing "x" stems...
#
    sumSIC = as.data.frame(table(sapply(lcs$stemsInCells, length)))
    colnames(sumSIC) = c('numStems', 'numCells')
    if(is.factor(sumSIC$numStems))
      sumSIC$numStems = as.integer(sumSIC$numStems)
      
#
#   get the data frame of the number of stems containing "x" cells...
#
    sumCIS = as.data.frame(table(sapply(lcs$cellsInStems, length)))
    colnames(sumCIS) = c('numCells', 'numStems')
    if(is.factor(sumCIS$numCells))
      sumCIS$numCells = as.integer(sumCIS$numCells)
      
#
#   get the cell and stem IDs just to make it easier for a user...
#
    cellIDs = names(lcs$stemsInCells)
    cellNums = as.integer(cellIDs)
    stemIDs = names(lcs$cellsInStems)
      
#
#   now create the object and return it...
#
    ssCSL = new('ssCellStemList',
                description = description,
                cellsInStems = lcs$cellsInStems,
                stemsInCells = lcs$stemsInCells,
                sumSIC = sumSIC,
                sumCIS = sumCIS,
                cellIDs = cellIDs,
                cellNums = cellNums,
                stemIDs = stemIDs
               )

    return(ssCSL)
}   #ssCellStemList
)   #setMethod


