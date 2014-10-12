# R script for plotting experimental results
require(tikzDevice)

# Set fonts to match Latex
  library(Cairo)
  mainfont <- "Garamond"
  CairoFonts(regular = paste(mainfont,"style=Regular",sep=":"),
             bold = paste(mainfont,"style=Bold",sep=":"),
             italic = paste(mainfont,"style=Italic",sep=":"),
             bolditalic = paste(mainfont,"style=Bold Italic,BoldItalic",sep=":"))
  pdf <- CairoPDF
  png <- CairoPNG
#

# R --vanilla -f plot-results.R --args <links.csv> <plain.csv>
#options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) {
  message("Usage: R --vanilla -f plot-results.R --args <perf.links.csv> <perf.plain.csv> <stoch.links.csv> <stoch.plain.csv>")
  quit()
}

perf_links_file=args[1]
perf_plain_file=args[2]
stoch_links_file=args[3]
stoch_plain_file=args[4]

set_margin <- function() {
  par(oma=c(0,0,0,0))
  par(mar=c(5,4,1.5,1))
}

data = list()
data[[1]]=read.csv(file=perf_links_file, as.is=T,row.names=1)
data[[2]]=read.csv(file=perf_plain_file, as.is=T,row.names=1)
data[[3]]=read.csv(file=stoch_links_file,as.is=T,row.names=1)
data[[4]]=read.csv(file=stoch_plain_file,as.is=T,row.names=1)

kmers=as.numeric(data[[1]]['kmer',])
cols=c('firebrick','dodgerblue', 'darkolivegreen', 'darkorchid')

fields=c('median','mean','N50','length','contigs','med_walk')
descriptions=c("Median contig length","Mean contig length","Contig N50",
               "Assembled length","Number of contigs", "Median walking distance")
path='plots/'

for(i in 1:length(fields)) {
  field = fields[i]
  description = descriptions[i]

  pdf(file=paste(path,field,'.pdf',sep=''),width=7,height=7)
  # tikz(paste(path,field,'.tex',sep=''),width=2.5,height=2.5)
  set_margin()

  m = 0;
  for(j in 1:4) { m = max(m, as.numeric(data[[j]][field,])); }

  plot(kmers,as.numeric(data[[1]][field,]),type='b',col=cols[1],
       ylim=c(0,m),axes=F,ylab=description,xlab="kmer",
       )#main=paste(description,"vs kmer-size"))
  axis(side = 2)
  axis(side = 1,at=kmers)

  for(j in 2:4) {
    points(kmers,as.numeric(data[[j]][field,]),type='b',col=cols[j])
  }

  s0=0; sn=0;
  for(j in 1:4) {
    s0 = s0 + as.numeric(data[[j]][field,1]);
    sn = sn + as.numeric(data[[j]][field,length(kmers)-1]);
  }
  lpos='topright';
  if(s0 < sn) { lpos='bottomright' }
  legend(lpos,c('Perfect Links','Perfect Plain','Stoch. Links','Stoch. Plain'),fill=cols)
  dev.off()
}
