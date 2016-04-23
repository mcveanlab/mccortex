#!/usr/bin/env Rscript --vanilla

#
# Dot plot of ref allele length vs sample allele length
#

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) {
  stop("Usage: Rscript --vanilla large-events-plot.R <stats.txt>\n")
}

file <- args[1]
#file <- 'bubbles50K/stats.txt'

title <- expression(paste(italic('klebsiella pneumoniae'),' large event sizes'))
xlabel <- 'Reference length (kbp)'
ylabel <- 'Sample length (kbp)'

r <- read.table(file,sep='\t',head=F,comment.char='#')
r <- r / 1000

# Get maximum and round to nearest 5
lim <- max(r[,1],r[,2])
lim <- floor((ceiling(lim)+4)/5)*5

pdf(file='kleb_large_events_R.pdf', width=6, height=6)
plot(r, xlab=xlabel, ylab=ylabel, main=title, xlim=c(0,lim), ylim=c(0,lim))
dev.off()

# With ggplot
library('ggplot2')
library('reshape')
library('scales')
library('plyr')

df <- data.frame(ref=r[,1], sample=r[,2])

p <- ggplot(df, aes(x=ref, y=sample)) +
     geom_point(shape=1) +
     xlim(0,lim) + ylim(0,lim) +
     ggtitle(title) + xlab(xlabel) + ylab(ylabel)

ggsave(p, file='kleb_large_events_ggplot.pdf', width=6, height=6)
