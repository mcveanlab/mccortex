#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2016-12-16

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) {
  stop("Usage: ./plot-mccortex-vs-sga.R <ng50.pdf> <errs.pdf> <mccortex.csv> <sga.csv>\n")
}


plot_path <- "latest/links-vs-sga-ng50.pdf"
errs_path <- "latest/links-vs-sga-errs.pdf"
mccortex_csv <- "latest/stocherr.links.csv"
sga_csv <- "latest/stocherr.sga.csv"

plot_path = args[1]
errs_path = args[2]
mccortex_csv <- args[3]
sga_csv <- args[4]

a <- read.table(mccortex_csv, sep=',',head=T,comment.char='#',as.is=T)
b <- read.table(sga_csv, sep=',',head=T,comment.char='#',as.is=T)

a$NG50 <- a$NG50/1000
b$NG50 <- b$NG50/1000

NG50_ylim <- ceiling(max(a$NG50,b$NG50)/20)*20
errs_ylim <- ceiling(max(b$AssemblyErrors)/100)*100
kmers <- a$K

# Plotting parameters
cols <- c('#1b9e77', '#d95f02', '#7570b3', 'red') # from color brewer
pnts <- c(19,4,17,1) # point styles pch=
jf <- 0.2 # jitter factor
lt <- 2.5 # line thickness
#

xlabel = expression(tau['min']~'(SGA) or'~italic('k')~'(McCortex)')

# pdf(plot_path, width=6, height=6)
quartz(type='pdf',file=plot_path,width=6,height=5)

# Remove empty title space
par(mar=c(4,4,4,2)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

plot(1, type='n', bty="n", axes=F,
     xlab=xlabel, ylab='NG50 (Kbp)',
     xlim=c(20,100), ylim=c(0,NG50_ylim))
points(jitter(a$K,jf), a$NG50, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=1)
points(jitter(b$K,jf), b$NG50, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=1)
axis(1, at=kmers)
axis(2, at=seq(0,NG50_ylim,20), las=2)

par(xpd=TRUE)
legend("topleft", bty="n", inset=c(0.2,-0.15),
       legend=c("McCortex","SGA"),
       col=cols[3:4], lwd=lt, lty=c(1,1,1),
       pch=pnts[3:4])

dev.off()



# pdf(errs_path, width=6, height=6)
quartz(type='pdf',file=errs_path,width=6,height=4)

# Remove empty title space
par(mar=c(4,4,4,2)+0.1) # set margins: bottom, left, top and right
par(xpd=TRUE)

# Concatenate with ~ (adds a space separator)
plot(1, type='n', bty="n", log="y", axes=F,
     xlab=xlabel, ylab='# Assembly errors (log)',
     xlim=c(20,100), ylim=c(1,errs_ylim))
points(jitter(a$K,jf), a$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[3], col=cols[3], lty=1)
points(jitter(b$K,jf), b$AssemblyErrors+1, type='b', lwd=lt, pch=pnts[4], col=cols[4], lty=1)
axis(1, at=kmers)
axis(2, at=c(1,2,11,51,101,201,501), labels=c(0,1,10,50,100,200,500), las=2)

par(xpd=TRUE)
legend("topright", bty="n", inset=c(0.2,-0.15),
       legend=c("McCortex","SGA"),
       col=cols[3:4], lwd=lt, lty=c(1,1,1),
       pch=pnts[3:4])

dev.off()
