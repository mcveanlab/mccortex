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

title <- expression(paste(italic('Klebsiella pneumoniae'),' large event sizes'))
xlabel <- 'Reference length (kbp) (log)'
ylabel <- 'Sample length (kbp) (log)'

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

lim <- max(lim_x, lim_y)
lim_x <- lim
lim_y <- lim

lim_x <- 100
lim_y <- 100
axis_labels = c(0.01, 0.1, 1, 10, 100)

pdf(file='kleb_large_events_combined_R_log.pdf', width=6, height=6)
plot(r1, xlab=xlabel, ylab=ylabel, log="xy", axes=FALSE,
     main=title, xlim=c(0.01,lim_x), ylim=c(0.01,lim_y))
points(r2)
axis(1, at=axis_labels, labels=axis_labels)
axis(2, at=axis_labels, labels=axis_labels)
dev.off()

# With ggplot
library('ggplot2')
library('reshape')
library('scales')
library('plyr')

shared <- data.frame(ref=r1[,1], sample=r1[,2])
unique <- data.frame(ref=r2[,1], sample=r2[,2])
shared['src'] <- 'shared'
unique['src'] <- 'links'
df <- rbind(shared, unique)

shared_label <- 'With or without links'
unique_label <- 'Links only'
legend_title <- "Events found by call set"

# Ideally shapes 20 and 4
p <- ggplot(df, aes(x=ref, y=sample, color=src, shape=src)) +
     geom_point(size=2) +
     scale_shape_manual(values=c('links'=17,'shared'=20),
                        labels=c('links'=unique_label, 'shared'=shared_label),
                        name=legend_title) +
     scale_color_manual(values=c('links'="#F8766D",'shared'="#00BFC4"),
                        labels=c('links'=unique_label, 'shared'=shared_label),
                        name=legend_title) +
     # scale_colour_brewer(name=legend_title) +
     # scale_color_manual(values=c("red", "blue"), name=legend_title) +
     scale_x_log10(limits=c(0.005, 150), breaks=axis_labels, labels=axis_labels) +
     scale_y_log10(limits=c(0.005, 150), breaks=axis_labels, labels=axis_labels) +
     # xlim(0, 60) + ylim(0, 60) +
     ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
     # scale_colour_discrete(name=legend_title) +
     # scale_shape_discrete(name=legend_title) +
     theme(plot.title = element_text(hjust = 0.5),
           legend.justification = c(1, 1), legend.position = c(1, 1),
           legend.background = element_rect(fill=alpha('white', 0.4)))

ggsave(p, file='kleb_large_events_combined_ggplot_log.pdf', width=6, height=6)
