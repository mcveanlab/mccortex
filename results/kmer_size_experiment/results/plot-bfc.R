#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2016-10-12

library('ggplot2')
library('reshape')
library('scales')
library('plyr')
library('gridExtra')
library('cowplot')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 8) {
  stop("Usage: ./plot-n50-and-errs.R <title> <out.pdf> <plain.csv> <links.csv> <pe.csv> <p.c.csv> <l.c.csv> <pe.c.csv>\n")
}

#
plain_csv <- "latest/stocherr.plain.csv"
links_csv <- "latest/stocherr.links.csv"
pe_csv <- "latest/stocherr.pe.csv"
plot_title <- "Stochastic coverage + 0.5% error (100X, 100bp reads)"
corr_plain_csv <- "latest/stocherrcorr.plain.csv"
corr_links_csv <- "latest/stocherrcorr.links.csv"
corr_pe_csv <- "latest/stocherrcorr.pe.csv"
output_pdf <- "plot.pdf"
#

use_pe <- (length(args) >= 5)
use_corr <- (length(args) >= 8)

plot_title <- args[1]
output_pdf <- args[2]
plain_csv <- args[3]
links_csv <- args[4]
if(use_pe) {
  pe_csv <- args[5]
}
if(use_corr) {
  corr_plain_csv <- args[6]
  corr_links_csv <- args[7]
  corr_pe_csv <- args[8]
}

a <- read.table(plain_csv,sep=',',head=T,comment.char='#',as.is=T)
b <- read.table(links_csv,sep=',',head=T,comment.char='#',as.is=T)
c <- read.table(pe_csv,sep=',',head=T,comment.char='#',as.is=T)
a$graph = factor('plain')
b$graph = factor('links')
c$graph <- factor('pe')

aa <- read.table(corr_plain_csv,sep=',',head=T,comment.char='#',as.is=T)
bb <- read.table(corr_links_csv,sep=',',head=T,comment.char='#',as.is=T)
cc <- read.table(corr_pe_csv,sep=',',head=T,comment.char='#',as.is=T)

a$corrNG50 <- aa$NG50
b$corrNG50 <- bb$NG50
c$corrNG50 <- cc$NG50
a$corrAsmErrors <- aa$AssemblyErrors
b$corrAsmErrors <- bb$AssemblyErrors
c$corrAsmErrors <- cc$AssemblyErrors

d = rbind(a,b,c)
dlevels <- c('pe','links','plain')
dlabels <- c('links PE','links','plain')
d$graph <- factor(d$graph, levels=dlevels, labels=dlabels)

# Approach 1
# Plot contig N50
p1 <- ggplot(data=d, aes(x=K, y=NG50, color=graph, shape=graph)) + theme_minimal() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
      geom_line() + geom_point(shape=4) +
      geom_line(aes(y=corrNG50),linetype="dotted") +
      scale_y_continuous(limits = c(0,250000)) +
      ylab("NG50") + ggtitle(plot_title) +
      theme(legend.title=element_blank()) + # hide legend title
      theme(legend.justification=c(0,1), legend.position=c(0,1)) # legend in plot top left

# Plot assembly error rate
p2 <- ggplot(data=d, aes(x=K, y=AssemblyErrors, color=graph)) + theme_minimal() +
      geom_point(shape=4) + geom_line() +
      geom_line(aes(y=corrAsmErrors),linetype="dotted") +
      scale_y_continuous(breaks=seq(0,150,50)) + coord_cartesian(ylim=c(0,150)) +
      ylab("Assembly Errors") +
      theme(legend.position="none") # hide legend

g <- plot_grid(p1, p2, align="v", nrow=2, rel_heights=c(3, 1))

ggsave(g, file=output_pdf, width=6, height=6)
