#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2017-02-16

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) {
  stop("Usage: ./plot-link-counts.R <linkcounts.pdf> <linkcounts.csv>\n")
}

plot_path <- "latest/perfect.linkcounts.se.pdf"
csv_path <- "latest/perfect.linkcounts.se.csv"

plot_path = args[1]
csv_path = args[2]

a <- read.table(csv_path, sep='\t',head=T,comment.char='#',as.is=T)

# Plotting parameters
cols <- c('#1b9e77', '#d95f02', '#7570b3', 'red') # from color brewer
pnts <- c(19,4,17,1) # point styles pch=
jf <- 0.2 # jitter factor
lt <- 2.5 # line thickness
#

# * joins with no spaces, ~ joins with a space
xlabel = expression(italic('k'))
ylabel = expression('no. of '*italic('k')*'mers with links (log)')

# pdf(plot_path, width=6, height=6)
quartz(type='pdf',file=plot_path,width=6,height=5)

# Remove empty title space
par(mar=c(4,5,2,2)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

par(mgp=c(4, 1, 0)) # axis label positions

plot(a$K, a$n_link_kmers, type='b', axes=F, log='y',
     xlab='', ylab='', ylim=c(1,max(a$n_link_kmers)))

mtext(side=1, text=xlabel, line=2)
mtext(side=2, text=ylabel, line=4)
axis(1, at=a$K)
axis(2, las=2)

dev.off()
