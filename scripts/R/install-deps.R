#!/usr/bin/env Rscript --vanilla

#
# Install all R packages required by McCortex R scripts
#

getpkg <- function(pkg) {
  if(!require(pkg, character.only=TRUE)) {
    install.packages(pkg, dep=TRUE, repos='http://cran.rstudio.com/')
  }
}

getpkg('ggplot2')
getpkg('gridExtra')
getpkg('reshape')
getpkg('scales')
getpkg('plyr')
