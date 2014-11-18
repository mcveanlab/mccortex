#script due to Richard Pearson and Zam

## arguments are filename (.covg file), expected depth of coverage,
## and a boolean (0/1) saying whether you want to print a pdf,
## and a boolean (0/1) saying whether you want to be in debug mode
## (will print graphs for D1 and D2 also)


args <- commandArgs(trailingOnly = TRUE) ##  one argument only -  start number
filename<-args[1]
depth<-as.integer(args[2])
printpdf<-as.integer(args[3])
debug<-as.integer(args[4])

threshold=1
print(filename)
dat <- read.table(filename, header=T)

##This is a bit like the first derivative (but is not)
d1 <- c(0, dat[,2]) / c(dat[,2], 0)
d1 <- d1[-c(1:2)]
#this is a bit like the 2nd derivative
d2 <- c(0, d1) / c(d1, 0)
d2 <- d2[-1]

if(printpdf==1)
{
  pdf(paste(filename, "pdf", sep="."))
  if(debug==1)
  {
    par(mfrow=c(1,3))
    plot(d1[1:(depth*2)], main="D1=Ratio of successive covgs", xlab="Kmer covg")
    plot(d2[1:(depth*2)], main="D2", xlab="Kmer covg")
  }
}

firstPosD1 <- min(which(d1 < 1))
firstPosD2 <- min(which(d2 <= threshold))

if(!is.na(firstPosD1) && (firstPosD1 < (depth * 0.75))) {
  cutoff <- firstPosD1
  #print("D1 cutoff:")
} else if(!is.na(firstPosD2)) {
  cutoff <- firstPosD2
  #print("D2 cutoff:")
} else {
  cutoff <- max(1, ceiling(depth/2))
  #print("Coverage cutoff:")
}

if(printpdf==1)
{
  plot(dat[, 1], dat[,2], log="y", xlim=c(1,depth*2),
       main=("Kmer coverage (with chosen cleaning threshold shown)"),
       xlab="Kmer covg", ylab="Frequency (log scale)")
  abline(v=cutoff)
  dev.off()
}

print("Returning cutoff is ")
print(cutoff)

