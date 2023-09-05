######################################################################
#
# InitErgmTerm.components.R
#
# copyright (c) 2022, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/02/23
# Licensed under the GNU General Public License v3 or later
#
# Portions based on the statnet ergm.userterms package (Hunter et al.,
# 2013, JSS).
#
# Part of the R/ergm.components package
#
# This file contains InitErgmTerm functions for component and
# connectivity related terms.
#
######################################################################


#InitErgmTerm function for the "bridges" term
#
#This term counts the number of bridging edges (cutedges) in the graph; in the
#directed case, weak connectivity is used.
InitErgmTerm.bridges <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
      varnames = "dyadic",
      vartypes = "logical",
      required = FALSE,
      defaultvalues = list(TRUE))
  if(is.directed(nw)){
    if(a$dyadic)
      cn<-"bridges.dyadic"
    else
      cn<-"bridges.cutedges"
  }else
    cn<-"bridges"
  list(name = "bridges",
      coef.names = cn,
      pkgname = "ergm.components",
      dependence = TRUE,
      emptynwstats = 0,
      inputs = a$dyadic,
      minval=0,
      maxval=network.size(nw)-1
  )
}


#InitErgmTerm function for the "components" term
#
#This term counts the number of components in the graph; in the directed case,
#these are taken to be weak components.
InitErgmTerm.components <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, 
      varnames = NULL,
      vartypes = NULL,
      required = NULL,
      defaultvalues = list())
  list(name = "components",
      coef.names = "components",
      pkgname = "ergm.components",
      dependence = TRUE,
      emptynwstats = network.size(nw),
      minval=1,
      maxval=network.size(nw)
  )
}


#InitErgmTerm function for the "compsizesum" term
#
#This term computes a function of the form
#
#  sum_c f(|V(c)|)
#
#where the sum is over the components (in the directed case, weak components)
#of the graph.  Currently, f is taken to be either
#
#  f(x) = x^k
#
#if log=FALSE, or
#
#  f(x) = log(x)^k
#
#if log=TRUE.  k is set by the "power" argument.
InitErgmTerm.compsizesum <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
      varnames = c("log", "power"),
      vartypes = c("logical","numeric"),
      required = c(FALSE,FALSE),
      defaultvalues = list(FALSE, 2))
  uselog<-a$log
  pow<-a$power
  cnam<-"compsizesum"
  if(uselog)
    cnam<-paste(cnam,"log",sep=".")
  if(uselog){
    minval<-0
    maxval<-NULL  #Unfortunately, not trivial in the general case
    if(pow==0)
      enws<-network.size(nw)
    else
      enws<-0
  }else{
    if(pow<0){
      minval<-network.size(nw)^pow
      maxval<-network.size(nw)
    }else{
      minval<-network.size(nw)
      maxval<-network.size(nw)^pow
    }
    enws<-network.size(nw)
  }
  cnam<-paste(cnam,"pow",pow,sep=".")
  list(name = "compsizesum",
      coef.names = cnam,
      pkgname = "ergm.components",
      dependence = TRUE,
      inputs = c(uselog,pow),
      emptynwstats = enws,
      minval=minval,
      maxval=maxval
  )
}


#InitErgmTerm function for the "dimers" term
#
#This term counts the number of components of order 2 in the graph; in the
#directed case, these are taken to be weak components.
InitErgmTerm.dimers <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
      varnames = NULL,
      vartypes = NULL,
      required = NULL,
      defaultvalues = list())
  if(is.directed(nw))
    fnam<-"dimers_dir"
  else
    fnam<-"dimers_udir"
  list(name = fnam,
      coef.names = "dimers",
      pkgname = "ergm.components",
      dependence = TRUE,
      emptynwstats = 0,
      minval=0,
      maxval=floor(network.size(nw)/2)
  )
}


