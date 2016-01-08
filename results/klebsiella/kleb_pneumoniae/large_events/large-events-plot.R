#!/usr/bin/env Rscript --vanilla

#
# Dot plot of ref allele length vs sample allele length
#

file='bubbles50K/stats.txt'

title <- expression(paste(italic('klebsiella pneumoniae'),' large event sizes'))
xlabel <- 'Reference length (kbp)'
ylabel <- 'Sample length (kbp)'

r <- read.table(file,sep='\t',head=F,comment.char='#')
r <- r / 1000

pdf(file='kleb_large_events_R.pdf', width=6, height=6)
plot(r, xlab=xlabel, ylab=ylabel, main=title, xlim=c(0,35), ylim=c(0,35))
dev.off()

# With ggplot
library('ggplot2')
library('reshape')
library('scales')
library('plyr')

df <- data.frame(ref=r[,1], sample=r[,2])

p <- ggplot(df, aes(x=ref, y=sample)) +
     geom_point(shape=1) +
     xlim(0,35) + ylim(0,35) +
     ggtitle(title) + xlab(xlabel) + ylab(ylabel)

ggsave(p, file='kleb_large_events_ggplot.pdf', width=6, height=6)
