#$Author: sinnwell $
#$Date: 2004/03/01 20:57:36 $
#$Header: /people/biostat3/sinnwell/Haplo/Make/RCS/setupData.q,v 1.1 2004/03/01 20:57:36 sinnwell Exp $
#$Locker:  $
#$Log: setupData.q,v $
#Revision 1.1  2004/03/01 20:57:36  sinnwell
#Initial revision
#

# use this function to define an alias function to run exactly as data() in R
# does nothing in Splus.
# R keeps a data set within the working data frame, so we only want to load
#     it when calling an example 
# Splus keeps it in background, so we have it loaded upon library(mypkg)

#=======================================
# Sinnwell JP
# Mayo Clinic, Rochester MN
# 2/28/2004
#=======================================

setupData <- function(...) {
  result <- if(exists("is.R") && is.function(is.R) && is.R()) {
               data
             } else {
               function(...) invisible(NULL)
             }
  eval(result(...))
}
