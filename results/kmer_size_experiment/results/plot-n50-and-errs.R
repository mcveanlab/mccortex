#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2016-10-12

library('ggplot2')
library('reshape')
library('scales')
library('plyr')
library('gridExtra')
library('cowplot')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) {
  stop("Usage: ./plot-n50-and-errs.R <plain.csv> <links.csv> <title> <out.pdf>\n")
}

# plain_csv <- "perfect.plain.csv"
# links_csv <- "perfect.links.csv"
# plot_title <- "Perfect coverage (100X, 100bp reads)"

# plain_csv <- "stoch.plain.csv"
# links_csv <- "stoch.links.csv"
# plot_title <- "Stochastic coverage (100X, 100bp reads)"

# plain_csv <- "stocher.plain.csv"
# links_csv <- "stocher.links.csv"
# plot_title <- "Stochastic coverage + Error (100X, 100bp reads, 1% err)"

# output_pdf <- "plot.pdf"

plain_csv <- args[1]
links_csv <- args[2]
plot_title <- args[3]
output_pdf <- args[4]

p <- read.table(plain_csv,sep=',',head=T,comment.char='#',as.is=T)
p$graph = 'plain'
l <- read.table(links_csv,sep=',',head=T,comment.char='#',as.is=T)
l$graph = 'links'
d <- rbind(p,l)

# Approach 1
p1 <- ggplot(data=d, aes(x=K, y=NG50, color=graph)) + theme_minimal() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
      geom_point(shape=4) + geom_line() +
      scale_y_continuous(limits = c(0,80000)) +
      ylab("NG50") + ggtitle(plot_title) +
      theme(legend.title=element_blank()) + # hide legend title
      theme(legend.justification=c(1,1), legend.position=c(1,1)) # legend in plot

p2 <- ggplot(data=d, aes(x=K, y=AssemblyErrors, color=graph)) + theme_minimal() +
      geom_point(shape=4) + geom_line() +
      scale_y_continuous(breaks=seq(0,20,10)) + coord_cartesian(ylim=c(0,20)) +
      ylab("Assembly Errors") +
      theme(legend.position="none") # hide legend

# 1a.
# grid.arrange(p1, p2, ncol=1, heights=c(2, 1))

# 1b.
g <- plot_grid(p1, p2, align="v", nrow=2, rel_heights=c(3, 1))

# 1c.
# grid.newpage()
# grid.draw(rbind(ggplotGrob(q), ggplotGrob(s), size = "last"))

# Approach 2
# m <- melt(d,measure.vars=c('NG50','AssemblyErrors'))

# ggplot(m, aes(x=K, y=value)) +
#   geom_line(aes(color=graph)) +
#   facet_grid(variable ~ ., scales="free_y", space="fixed") +
#   xlab("kmer size") + ylab("") + ggtitle(plot_title)

ggsave(g, file=output_pdf, width=6, height=6)
