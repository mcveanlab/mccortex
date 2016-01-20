#!/usr/bin/env Rscript --vanilla


load_data <- function(file) {
  r <- read.table(file,sep=',',head=F,comment.char='#',as.is=T)
  return(r)
}

hists <- list()
hists[[1]] <- load_data('hists/platypus.hist.csv')
hists[[2]] <- load_data('hists/freebayes.hist.csv')
hists[[3]] <- load_data('hists/cortex.hist.csv')
hists[[4]] <- load_data('hists/bubbles.hist.csv')
hists[[5]] <- load_data('hists/breakpoints.hist.csv')

names <- c('Platypus','Freebayes','Cortex', 'McCortex Bubbles', 'McCortex Breakpoints')
colours <- c('red', 'green', 'blue', 'black', 'purple')
# points <- c('o','o','o', 'o', 'o')
points <- rep(20, 5)

title <- expression(paste(italic('klebsiella pneumoniae'),' indels'))

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

# plot

jf=2
pdf(file="kleb_indels.pdf",width=6,height=6)

plot(jitter(x,factor=jf), m[,2], log="y",
     col=colours[1],pch=points[1],
     main=title, xlab="indel size", ylab="count (log)")

for(i in 2:length(names)) {
  points(jitter(x,factor=jf),m[,i+1],col=colours[i],pch=points[i]) # ,type='b'
}
legend("topright", names, col=colours, pch=points)

dev.off()

# dm <- data.frame(m)
# p <- ggplot(dm, aes(x=dist, y=Platypus)) + geom_line()
