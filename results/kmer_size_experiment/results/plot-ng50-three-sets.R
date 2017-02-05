#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2016-12-16

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 7) {
  stop("Usage: ./plot-n50-and-errs.R <out.pdf> [<X.plain.csv> <X.links.csv>]x['perfect','stoch','stocherr']\n")
}


plot_path = "plot.pdf"
perf_plain_csv <- "latest/perfect.plain.csv"
perf_links_csv <- "latest/perfect.pe.csv"
stoch_plain_csv <- "latest/stoch.plain.csv"
stoch_links_csv <- "latest/stoch.pe.csv"
error_plain_csv <- "latest/stocherr.plain.csv"
error_links_csv <- "latest/stocherr.pe.csv"

plot_path = args[1]
perf_plain_csv <- args[2]
perf_links_csv <- args[3]
stoch_plain_csv <- args[4]
stoch_links_csv <- args[5]
error_plain_csv <- args[6]
error_links_csv <- args[7]

a <- read.table(perf_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
b <- read.table(perf_links_csv, sep=',',head=T,comment.char='#',as.is=T)
c <- read.table(stoch_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
d <- read.table(stoch_links_csv, sep=',',head=T,comment.char='#',as.is=T)
e <- read.table(error_plain_csv, sep=',',head=T,comment.char='#',as.is=T)
f <- read.table(error_links_csv, sep=',',head=T,comment.char='#',as.is=T)

a$NG50 <- a$NG50/1000
b$NG50 <- b$NG50/1000
c$NG50 <- c$NG50/1000
d$NG50 <- d$NG50/1000
e$NG50 <- e$NG50/1000
f$NG50 <- f$NG50/1000

NG50_ylim <- ceiling(max(a$NG50,b$NG50,c$NG50,d$NG50,e$NG50,f$NG50)/20)*20
kmers <- a$K

xlabel = expression(italic('k'))

pdf(plot_path, width=6, height=6)

# Remove empty title space
par(mar=c(4,4,4,2)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

cols=c('#1b9e77', '#d95f02', '#7570b3') # from color brewer
pnts <- c(19,4,17) # point styles pch=
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
axis(1, at=kmers)
axis(2, at=seq(0,NG50_ylim,20), las=2)

par(xpd=TRUE)
legend("topleft", bty="n", inset=c(0.2,-0.15),
       legend=c("perfect","stochastic","error"),
       col=cols, lwd=lt, lty=c(1,1,1), 
       pch=pnts)

legend("topright", -0.2, bty="n", inset=c(0.2,-0.15),
       legend=c("with links","without links"),
       col="black", lwd=lt, lty=c(1,2), pch='')

dev.off()


# ggplot2 version

# a$src <- factor('perfect')
# b$src <- factor('perfect')
# c$src <- factor('stoch')
# d$src <- factor('stoch')
# e$src <- factor('stocherr')
# f$src <- factor('stocherr')

# a$graph <- factor('plain')
# b$graph <- factor('links')
# c$graph <- factor('plain')
# d$graph <- factor('links')
# e$graph <- factor('plain')
# f$graph <- factor('links')

# g <- rbind(a,b,c,d,e,f)

# p1 <- ggplot(data=g, aes(x=K, y=NG50)) +
#       # geom_line(aes(group=interaction(src, graph), colour=src, linetype=graph), linejoin='mitre') +
#       geom_point(aes(shape=src, color=src)) +
#       theme(legend.title=element_blank()) # hide legend title

#       geom_line(aes(linetype='dotted'))
#       geom_line(graph, aes(color=src, linetype=graph)) + geom_point(aes(shape=src)) +
#       scale_y_continuous(limits = c(0,NG50_ylim)) +
#       ylab("NG50 (kbp)") + ggtitle(plot_title) +
#       theme(legend.title=element_blank()) + # hide legend title
#       theme(legend.justification=c(0,1), legend.position=c(0,1)) # legend in plot top left

