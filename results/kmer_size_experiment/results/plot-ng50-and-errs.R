#!/usr/bin/env Rscript --vanilla

# Isaac Turner 2016-10-12

library('ggplot2')
library('reshape')
library('scales')
library('plyr')
library('gridExtra')
library('cowplot')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4 && length(args) != 5 && length(args) != 8) {
  stop("Usage: ./plot-ng50-and-errs.R <title> <out.pdf> <plain.csv> <links.csv> [[pe.csv] [<p.c.csv> <l.c.csv> <pe.c.csv>]]\n")
}

# plain_csv <- "perfect.plain.csv"
# links_csv <- "perfect.links.csv"
# pe_csv <- "perfect.pe.csv"
# plot_title <- "Perfect coverage (100X, 100bp reads)"
# use_pe <- TRUE
# use_corr <- FALSE

# plain_csv <- "latest/stoch.plain.csv"
# links_csv <- "latest/stoch.links.csv"
# pe_csv <- "latest/stoch.pe.csv"
# plot_title <- "Stochastic coverage (100X, 100bp reads)"
# use_pe <- TRUE
# use_corr <- FALSE

plain_csv <- "latest/stocherr.plain.csv"
links_csv <- "latest/stocherr.links.csv"
pe_csv <- "latest/stocherr.pe.csv"
plot_title <- "Stochastic coverage + 0.5% err (100X, 100bp reads)"
use_pe <- TRUE
use_corr <- FALSE

plain_csv <- "latest/stocherr.plain.csv"
links_csv <- "latest/stocherr.links.csv"
pe_csv <- "latest/stocherr.pe.csv"
plot_title <- "Stochastic coverage + 0.5% error (100X, 100bp reads)"
corr_plain_csv <- "latest/stocherrcorr.plain.csv"
corr_links_csv <- "latest/stocherrcorr.links.csv"
corr_pe_csv <- "latest/stocherrcorr.pe.csv"
use_pe <- TRUE
use_corr <- TRUE

# output_pdf <- "plot.pdf"

#
# Parse params
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
a$graph = factor('plain')
b <- read.table(links_csv,sep=',',head=T,comment.char='#',as.is=T)
b$graph = factor('links')

if(use_pe) {
  c <- read.table(pe_csv,sep=',',head=T,comment.char='#',as.is=T)
  c$graph <- factor('pe')
  dlevels <- c('pe','links','plain')
  dlabels <- c('links PE','links','plain')
} else {
  dlevels <- c('links','plain')
  dlabels <- c('links','plain')
}

if(use_corr) {
  aa <- read.table(corr_plain_csv,sep=',',head=T,comment.char='#',as.is=T)
  bb <- read.table(corr_links_csv,sep=',',head=T,comment.char='#',as.is=T)
  cc <- read.table(corr_pe_csv,sep=',',head=T,comment.char='#',as.is=T)
  a$corrNG50 <- aa$NG50 / 1000
  b$corrNG50 <- bb$NG50 / 1000
  c$corrNG50 <- cc$NG50 / 1000
  a$corrAsmErrors <- aa$AssemblyErrors
  b$corrAsmErrors <- bb$AssemblyErrors
  c$corrAsmErrors <- cc$AssemblyErrors
}

if(use_pe) {
  d <- rbind(a,b,c)
} else {
  d <- rbind(a,b)
}

d$graph <- factor(d$graph, levels=dlevels, labels=dlabels)
d$NG50 <- d$NG50 / 1000

N50_ylim <- 100
asm_ylim <- 5
if(use_pe) {
  N50_ylim <- 150
  asm_ylim <- 1000
}

xlabel = expression(italic('k'))

# Plot contig N50
p1 <- ggplot(data=d, aes(x=K, y=NG50, color=graph)) + theme_minimal() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
      geom_line() + geom_point(shape=4) +
      scale_x_continuous(breaks=d$K, minor_breaks=NULL) +
      scale_y_continuous(limits = c(0,N50_ylim)) +
      ylab("NG50 (kbp)") + xlab(xlabel) + ggtitle(plot_title) +
      theme(legend.title=element_blank()) + # hide legend title
      theme(legend.justification=c(0,1), legend.position=c(0,1)) # legend in plot top left

if(use_corr) {
  p1 <- p1 + geom_line(aes(y=corrNG50),linetype="dotted")
}

# Plot assembly error rate
p2 <- ggplot(data=d, aes(x=K, y=AssemblyErrors, color=graph)) + theme_minimal() +
      geom_line() + geom_point(shape=4) +
      scale_x_continuous(breaks=d$K, minor_breaks=NULL) +
      scale_y_continuous(limits=c(0, asm_ylim), trans="log1p") +
      ylab("Assembly Errors (log)") +
      theme(legend.position="none") # hide legend

if(use_corr) {
  p2 <- p2 + geom_line(aes(y=corrAsmErrors),linetype="dotted")
}

g <- plot_grid(p1, p2, align="v", nrow=2, rel_heights=c(3, 1))

ggsave(g, file=output_pdf, width=6, height=6)
