#---------------------------------------------------------------------------
#
#   This file holds the S4 class definitions for the "ssCellStemList" 
#   class.
#
#   Classes and subclasses in this file include...
#
#     1. ssCellStemList
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




#=================================================================================================
#
#  1. define the ssCellStemList class...
#
setClass('ssCellStemList',
#
#  slots for the class and its subclasses...
#
    representation(description = 'character',
                   cellsInStems = 'list',       #the list of cells w/in stem IZs
                   stemsInCells = 'list',       #the list of stems w/in cells
                   sumSIC = 'data.frame',       #summary data frame: number of cells with x stems (IZs)
                   sumCIS = 'data.frame',       #summary data frame: number of stems (IZs) with x cells
                   cellIDs = 'character',       #character vector of the cell IDs in the SIC list names
                   cellNums = 'integer',        #integer vector of cellIDs
                   stemIDs = 'character'        #character vector of stem IDs from the CIS list names
                   ),
                   
    #contains = 'VIRTUAL',                                      #note!
    
    #some defaults for validity checking...
    prototype = list(description = 'list of cells in stems & stems in cells',
                     cellsInStems = list(),     #dummy
                     stemsInCells = list(),     #dummy
                     sumSIC = data.frame(),     #dummy
                     sumCIS = data.frame(),     #dummy
                     cellIDs = character(0),
                     cellNums = integer(0),
                     stemIDs = character(0)
                    ),
                    
    validity = function(object) {

#              no checks for now...
                   
                 return(TRUE)
               } #validity check
) #class ssCellStemList 



