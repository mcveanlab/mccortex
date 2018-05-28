#!/usr/bin/env Rscript --vanilla

#
# Dot plot of ref allele length vs sample allele length
#

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) {
  stop("Usage: Rscript --vanilla single-plot.R <shared.txt> <links.txt>\n")
}

file1 <- args[1]
file2 <- args[2]
#file <- 'bubbles50K/stats.txt'

title <- expression(paste(italic('klebsiella pneumoniae'),' large event sizes'))
xlabel <- 'Reference length (kbp)'
ylabel <- 'Sample length (kbp)'

# Read files (format: <ref_len><tab><alt_len>), convert to kbp
r1 <- read.table(file1,sep='\t',head=F,comment.char='#')
r1 <- r1 / 1000

r2 <- read.table(file2,sep='\t',head=F,comment.char='#')
r2 <- r2 / 1000

# Get maximum and round to nearest 5
# lim_x <- max(r1[,1],r2[,1],r1[,2],r2[,2])
# lim_x <- floor((ceiling(lim_x)+4)/5)*5
# lim_y <- lim_x

lim_x <- max(r1[,1],r2[,1])
lim_x <- floor((ceiling(lim_x)+4)/5)*5
lim_y <- max(r1[,2],r2[,2])
lim_y <- floor((ceiling(lim_y)+4)/5)*5

# pdf(file='kleb_large_events_combined_R.pdf', width=6, height=6)
# plot(r1, xlab=xlabel, ylab=ylabel, main=title, xlim=c(0,lim), ylim=c(0,lim))
# points(r2)
# dev.off()

# With ggplot
library('ggplot2')
library('reshape')
library('scales')
library('plyr')

shared <- data.frame(ref=r1[,1], sample=r1[,2])
unique <- data.frame(ref=r2[,1], sample=r2[,2])
shared['src'] <- 'shared'
unique['src'] <- 'links only'
df <- rbind(shared, unique)

legend_title <- "Call set"

# Ideally shapes 20 and 4
p <- ggplot(df, aes(x=ref, y=sample, shape=src, color=src)) +
     geom_point(size=2) +
     scale_shape_manual(values=c(17,1), name=legend_title) +
     xlim(0, lim_x) + ylim(0, lim_y) +
     ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
     scale_colour_discrete(name=legend_title) +
     # scale_shape_discrete(name=legend_title) +
     theme(plot.title = element_text(hjust = 0.5),
           legend.justification = c(1, 1), legend.position = c(1, 1))

ggsave(p, file='kleb_large_events_combined_ggplot.pdf', width=6, height=6)
