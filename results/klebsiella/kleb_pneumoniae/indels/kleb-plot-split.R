#!/usr/bin/env Rscript --vanilla


load_data <- function(file) {
  r <- read.table(file,sep=',',head=F,comment.char='#',as.is=T)
  return(r)
}

hists <- list()
hists[[1]] <- load_data('hists/bubbles.hist.csv')
hists[[2]] <- load_data('hists/breakpoints.hist.csv')
hists[[3]] <- load_data('hists/platypus.hist.csv')
hists[[4]] <- load_data('hists/freebayes.hist.csv')
hists[[5]] <- load_data('hists/cortex.hist.csv')

names <- c('McCortex Bubbles', 'McCortex Breakpoints', 'Platypus','Freebayes','Cortex')
# colours <- c('red', 'green', 'blue', 'black', 'purple')

library(RColorBrewer)
colours <- brewer.pal(5, "Set1")
# colours <- rep('black', 5)


title <- expression(paste(italic('Klebsiella pneumoniae'),' indels'))

ncalls <- length(hists)
lim <- 100
m <- matrix(0, nrow=2*lim+1, ncol=ncalls+1)
m[,1] <- -lim:lim
colnames(m) <- c("dist",names)
x <- -lim:lim

# Merge
for(i in 1:length(hists)) {
  l <- hists[[i]]
  for(j in 1:nrow(l)) {
    d <- l[j,1]
    d <- min(d, lim)
    d <- max(d,-lim)
    d <- d + lim + 1
    m[d,1+i] <- m[d,1+i] + l[j,2]
  }
}

pdf(file="kleb_indels_split.pdf", width=10, height=10)

# attempting to force equal sized plot dimensions
par(mfrow=c(length(names), 1), cex=0.8, cex.lab=0.8, cex.axis=0.8, cex.main=1,
    mgp=c(1.5, 0.3, 0), oma=c(2, 1, 2, 1), pin=c(6,2), mar=c(1, 3, 1, 3),
    tcl=-0.2)

for(i in 1:length(names)) {
  plot(x, m[,1+i], log="y", col=colours[i], type="h",
       main=NA, xlab=NA, ylab="Count (log)",
       xlim=c(-lim, lim), ylim=c(1, max(m)))
  legend("topright", names[i], col=colours[i], bty="n")
}

title(title, outer=TRUE)
mtext("Indel size (bp)", side=1, line=0.5, outer=TRUE, cex=0.8)

dev.off()
