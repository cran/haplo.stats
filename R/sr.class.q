## $[id]

## Author: Jason Sinnwell
## Date: 10/2005
## Purpose: part of the srlocal strategic grant
##          Mayo Division of Biostatistics
##          an SR-universal function that sets all S3 classes
##          for inheritence of x

## oldClass from splus:
#       the class of an object as defined in Version 3 of S;
# class from R
#       R objects have a 'class' attribute, a character vector giving
#       the names of the classes which the object "inherits" from.

"sr.class<-" <- function(x, value) {
# value is a character vector to assign S3 class
# inheritance vector for object x

  if(is.R()) {
    class(x) <- value
  } else {
    setOldClass(value)
    oldClass(x) <- value[1]
  }
  x
}

sr.class <- function(x) {
# return S3 class inheritance vector for object x

  if(is.R())
    class(x) 
  else
    oldClass(x) 
}
