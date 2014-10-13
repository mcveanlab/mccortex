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
if(length(args) != 6) {
  message("Usage: R --vanilla -f plot-results.R --args <{perf,stoch,stocherr}.{links,paths}.csv>")
  quit()
}

files=args

set_margin <- function() {
  par(oma=c(0,0,0,0))
  par(mar=c(5,4,1.5,1))
}

data = list()
for(i in 1:6) {
  data[[i]]=read.csv(file=files[i], as.is=T,row.names=1)
  cat("files[",i,"] = '",files[i],"'\n",sep='')
}

kmers=as.numeric(data[[1]]['kmer',])
cols=c('firebrick','dodgerblue', 'darkolivegreen', 'darkorchid', 'orange', 'black')

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
  for(j in 1:6) { m = max(m, as.numeric(data[[j]][field,])); }

  plot(kmers,as.numeric(data[[1]][field,]),type='b',col=cols[1],
       ylim=c(0,m),axes=F,ylab=description,xlab="kmer",
       )#main=paste(description,"vs kmer-size"))
  axis(side = 2)
  axis(side = 1,at=kmers)

  for(j in 2:6) {
    points(kmers,as.numeric(data[[j]][field,]),type='b',col=cols[j])
  }

  s0=0; sn=0;
  for(j in 1:6) {
    s0 = s0 + as.numeric(data[[j]][field,1]);
    sn = sn + as.numeric(data[[j]][field,length(kmers)-1]);
  }
  lpos='topright';
  if(s0 < sn) { lpos='bottomright' }
  legtxt = c('Perfect Links', 'Perfect Plain',
             'Stoch. Links', 'Stoch. Plain',
             'Stoch. Err Links', 'Stoch. Err Plain');
  legend(lpos,legtxt,fill=cols)
  dev.off()
}
