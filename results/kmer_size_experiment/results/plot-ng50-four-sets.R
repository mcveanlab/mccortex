#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2016-12-16

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 10) {
  stop("Usage: ./plot-n50-and-errs.R <ng50.pdf> <errs.pdf> [<X.plain.csv> <X.links.csv>]x['perfect','stoch','stocherr','corr']\n")
}


plot_ng50_path = "plot.pdf"
plot_err_path = "plot.pdf"
perf_plain_csv <- "latest/perfect.plain.csv"
perf_links_csv <- "latest/perfect.pe.csv"
stoch_plain_csv <- "latest/stoch.plain.csv"
stoch_links_csv <- "latest/stoch.pe.csv"
error_plain_csv <- "latest/stocherr.plain.csv"
error_links_csv <- "latest/stocherr.pe.csv"
corr_plain_csv <- "latest/stocherrcorr.plain.csv"
corr_links_csv <- "latest/stocherrcorr.pe.csv"

plot_ng50_path = args[1]
plot_err_path = args[2]
perf_plain_csv <- args[3]
perf_links_csv <- args[4]
stoch_plain_csv <- args[5]
stoch_links_csv <- args[6]
error_plain_csv <- args[7]
error_links_csv <- args[8]
corr_plain_csv <- args[9]
corr_links_csv <- args[10]

a <- read.table(perf_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
b <- read.table(perf_links_csv, sep=',',head=T,comment.char='#',as.is=T)
c <- read.table(stoch_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
d <- read.table(stoch_links_csv, sep=',',head=T,comment.char='#',as.is=T)
e <- read.table(error_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
f <- read.table(error_links_csv, sep=',',head=T,comment.char='#',as.is=T)
g <- read.table(corr_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
h <- read.table(corr_links_csv, sep=',',head=T,comment.char='#',as.is=T)

a$NG50 <- a$NG50/1000
b$NG50 <- b$NG50/1000
c$NG50 <- c$NG50/1000
d$NG50 <- d$NG50/1000
e$NG50 <- e$NG50/1000
f$NG50 <- f$NG50/1000
g$NG50 <- g$NG50/1000
h$NG50 <- h$NG50/1000

NG50_ylim <- ceiling(max(a$NG50, b$NG50,
                         c$NG50, d$NG50,
                         e$NG50, f$NG50,
                         g$NG50, h$NG50) / 20) * 20
errs_ylim <- ceiling(max(a$AssemblyErrors, b$AssemblyErrors,
                         c$AssemblyErrors, d$AssemblyErrors,
                         e$AssemblyErrors, f$AssemblyErrors,
                         g$AssemblyErrors, h$AssemblyErrors)/20)*20
kmers <- a$K

xlabel = expression(italic('k'))

pdf(plot_ng50_path, width=6, height=6)

# Remove empty title space
par(mar=c(4,4,3,0.5)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

cols=c('#1b9e77', '#d95f02', '#7570b3', 'darkgray') # from color brewer
pnts <- c(19, 4, 17, 2) # point styles pch=
jf <- 0.2 # jitter factor
lt <- 2.5 # line thickness

plot(1, type='n', bty="n", xlab=xlabel, ylab='NG50 (Kbp)',
     xlim=c(20,100), ylim=c(0,NG50_ylim), axes=F)
points(jitter(a$K,jf), a$NG50, type='b', lwd=lt, pch=pnts[1], col=cols[1], lty=2)
points(jitter(b$K,jf), b$NG50, type='b', lwd=lt, pch=pnts[1], col=cols[1], lty=1)
points(jitter(c$K,jf), c$NG50, type='b', lwd=lt, pch=pnts[2], col=cols[2], lty=2)
points(jitter(d$K,jf), d$NG50, type='b', lwd=lt, pch=pnts[2], col=cols[2], lty=1)
points(jitter(e$K,jf), e$NG50, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=2)
points(jitter(f$K,jf), f$NG50, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=1)
points(jitter(g$K,jf), g$NG50, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=2)
points(jitter(h$K,jf), h$NG50, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=1)
axis(1, at=kmers)
axis(2, at=seq(0,NG50_ylim,20), las=2)

par(xpd=TRUE)
legend("topleft", bty="n", inset=c(0.2,-0.15),
       legend=c("perfect", "stochastic", "error", "corrected"),
       col=cols, lwd=lt, lty=c(1,1,1,1), 
       pch=pnts)

legend("topright", -0.2, bty="n", inset=c(0.2,-0.15),
       legend=c("with links","without links"),
       col="black", lwd=lt, lty=c(1,2), pch='')

dev.off()


#
# Plot error plot
#

pdf(plot_err_path, width=6, height=4)

# Remove empty title space
par(mar=c(4,4,3,0.5)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

cols=c('#1b9e77', '#d95f02', '#7570b3', 'darkgray') # from color brewer
pnts <- c(19, 4, 17, 2) # point styles pch=
jf <- 0.2 # jitter factor
lt <- 2.5 # line thickness

plot(1, type='n', bty="n", xlab=xlabel, ylab='# Assembly errors (log)',
     xlim=c(20,100), ylim=c(1, errs_ylim), axes=F, log='y')
points(jitter(a$K,jf), a$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[1], col=cols[1], lty=2)
points(jitter(b$K,jf), b$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[1], col=cols[1], lty=1)
points(jitter(c$K,jf), c$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[2], col=cols[2], lty=2)
points(jitter(d$K,jf), d$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[2], col=cols[2], lty=1)
points(jitter(e$K,jf), e$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=2)
points(jitter(f$K,jf), f$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=1)
points(jitter(g$K,jf), g$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=2)
points(jitter(h$K,jf), h$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=1)
axis(1, at=kmers)
axis(2, at=c(1, 11, 51, 101, 201, 501, 1001),
     labels=c(0, 10, 50, 100, 200, 500, 1000), las=2)

par(xpd=TRUE)
legend("topleft", bty="n", inset=c(0.3,-0.15),
       legend=c("perfect", "stochastic", "error", "corrected"),
       col=cols, lwd=lt, lty=c(1,1,1,1), 
       pch=pnts)

legend("topright", -0.2, bty="n", inset=c(0.1,-0.15),
       legend=c("with links", "without links"),
       col="black", lwd=lt, lty=c(1,2), pch='')

dev.off()
