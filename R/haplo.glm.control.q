#$Author: sinnwell $
#$Date: 2004/03/02 16:34:07 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo.glm.control.q,v 1.3 2004/03/02 16:34:07 sinnwell Exp $
#$Locker:  $
#$Log: haplo.glm.control.q,v $
#Revision 1.3  2004/03/02 16:34:07  sinnwell
#change T to TRUE
#
#Revision 1.2  2003/12/08 19:37:28  sinnwell
# changed F,T to FALSE,TRUE
#
#Revision 1.1  2003/09/16 16:01:29  schaid
#Initial revision
#
haplo.glm.control <- function(haplo.effect="add",
                              haplo.base = NULL,
                              haplo.freq.min=0.001,
                              sum.rare.min=0.001,
                              haplo.min.info=0.001,
                              keep.rare.haplo=TRUE,
                              glm.c=glm.control(maxit=500),
                              em.c=haplo.em.control()){


  chk <- charmatch(haplo.effect, c("additive", "dominant", "recessive"))
  if(is.na(chk)) stop("Invalid haplo.effect")


  if(haplo.freq.min < 0 | haplo.freq.min > .9) {
    warning("The value of haplo.freq.min is out of range, the devault value of 0.001 is used instead")
    haplo.freq.min <- 0.001
  }


  if(sum.rare.min < 0 | sum.rare.min > .9) {
    warning("The value of sum.rare.min is out of range, the devault value of 0.001 is used instead")
    sum.rare.min <- 0.001
  }

 if(haplo.min.info < 0 | haplo.min.info > .9) {
    warning("The value of haplo.min.info is out of range, the devault value of 0.001 is used instead")
    haplo.min.info <- 0.001
  }


  if(keep.rare.haplo!=TRUE & keep.rare.haplo!=FALSE){
    warning("The value of keep.rare.haplo is invalid, the devalut value of TRUE is used instead")
    keep.rare.haplo=TRUE
  }


  return(list(haplo.effect=haplo.effect,
              haplo.base = haplo.base,
              haplo.freq.min=haplo.freq.min,
              sum.rare.min=sum.rare.min,
              haplo.min.info=haplo.min.info,
              keep.rare.haplo=keep.rare.haplo,
              epsilon=glm.c$epsilon,
              maxit=glm.c$maxit,
              trace=glm.c$trace,
              em.control=em.c))
}

