#!/usr/bin/env Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 1) {
  stop("Usage: stats.R <times.csv>");
}

cat('Reading:',args[1],'\n');

x=read.csv(file=args[1],as.is=T,header=F)

cat('Rows:',nrow(x),'\n');

i=1
j=1

while(i <= nrow(x)) {
  cat(j,' [',i,':',i+4,'] mean: ',mean(x[i:(i+4),1]),' stddev: ',sd(x[i:(i+4),1]),'\n',sep='');
  i=i+5; j=j+1;
}
