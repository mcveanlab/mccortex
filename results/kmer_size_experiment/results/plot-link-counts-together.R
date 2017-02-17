#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2017-02-16

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 5) {
  stop("Usage: ./plot-link-counts.R <linkcounts.pdf> <perfect.csv> <stoch.csv> <stocherr.csv> <stocherrcorr.csv>\n")
}


plot_path <- "latest/linkcounts.se.pdf"
perf_path <- "latest/perfect.linkcounts.se.csv"
stoch_path <- "latest/stoch.linkcounts.se.csv"
error_path <- "latest/stocherr.linkcounts.se.csv"
errcorr_path <- "latest/stocherrcorr.linkcounts.se.csv"

plot_path = args[1]
perf_path <- args[2]
stoch_path <- args[3]
error_path <- args[4]
errcorr_path <- args[5]

a <- read.table(perf_path, sep='\t',head=T,comment.char='#',as.is=T)
b <- read.table(stoch_path, sep='\t',head=T,comment.char='#',as.is=T)
c <- read.table(error_path, sep='\t',head=T,comment.char='#',as.is=T)
d <- read.table(errcorr_path, sep='\t',head=T,comment.char='#',as.is=T)

nlinkkmers_ylim <- ceiling(max(a$n_link_kmers, b$n_link_kmers,
                               c$n_link_kmers, d$n_link_kmers))
kmers <- a$K

# * joins with no spaces, ~ joins with a space
xlabel = expression(italic('k'))
ylabel = expression('no. of '*italic('k')*'mers with links (log)')

pdf(plot_path, width=6, height=6)

# Remove empty title space
par(mar=c(4,5,2,2)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

cols=c('#1b9e77', '#d95f02', '#7570b3', 'darkgray') # from color brewer
pnts <- c(19,4,17,2) # point styles pch=
jf <- 0.2 # jitter factor
lt <- 2 # line thickness

plot(1, type='n', bty="n", xlab='', ylab='', log='y',
     xlim=c(20,100), ylim=c(1,nlinkkmers_ylim), axes=F)

points(jitter(a$K,jf), a$n_link_kmers, type='b', lwd=lt, pch=pnts[1], col=cols[1], lty=1)
points(jitter(b$K,jf), b$n_link_kmers, type='b', lwd=lt, pch=pnts[2], col=cols[2], lty=1)
points(jitter(d$K,jf), d$n_link_kmers, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=1)
points(jitter(c$K,jf), c$n_link_kmers, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=1)

mtext(side=1, text=xlabel, line=2)
mtext(side=2, text=ylabel, line=4)
axis(1, at=a$K)
axis(2, las=2)

par(xpd=TRUE)
legend("topright", bty="n", inset=c(0.2,0),
       legend=c("perfect","stochastic","error","corrected"),
       col=cols, lwd=lt, lty=c(1,1,1,1), 
       pch=pnts)

dev.off()
