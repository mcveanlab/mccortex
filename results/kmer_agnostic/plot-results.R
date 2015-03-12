# R script for plotting experimental results
# require(tikzDevice)

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
# nfiles=9
args <- commandArgs(trailingOnly = TRUE)
# if(length(args) != nfiles) {
#   message("Usage: R --vanilla -f plot-results.R --args <{perf,stoch,stocherr}.{plain,links,string}.csv>")
#   quit()
# }

out_path=args[1]
files=args[2:length(args)]
nfiles=length(files)

cat('out_path=',out_path,'\n');

titles=c()
for(i in 1:nfiles) { titles[i] = substr(files[i],1,nchar(files[i])-10) }

set_margin <- function() {
  par(oma=c(0,0,0,0))
  par(mar=c(5,4,1.5,1))
}

data = list()
for(i in 1:nfiles) {
  data[[i]]=read.table(file=files[i], sep=',', as.is=T, row.names=1, header=T)
  cat("files[",i,"] = '",files[i],"'\n",sep='')
}

kmers=as.numeric(data[[1]]['kmer',])

# tmpcols=rainbow(3)
# cols=tmpcols[rep(1:3,each=3)]
# linetypes=rep(c(1,2,3),times=3)
cols=rainbow(nfiles)
linetypes=rep(c(1,2,3),length.out=nfiles)
pchtypes=rep(c(0,3,4),length.out=nfiles)

fields=c('median','mean','N50','sum_length','contigs','med_walk')
descriptions=c("Median contig length","Mean contig length","Contig N50",
               "Assembled length","Number of contigs", "Median walking distance")

# Correct length field to account for kmer overlap
for(f in 1:nfiles) {
  data[[f]]['sum_length',] = data[[f]]['length',] -
                             pmax(as.numeric(data[[f]]['contigs',])-1, 0) *
                             (kmers-1);
}

# Generate multi-page pdf
pdf(file=out_path,width=7,height=7)

for(i in 1:length(fields)) {
  cat('run',i,'\n');
  field = fields[i]
  description = descriptions[i]

  # tikz(paste(path,field,'.tex',sep=''),width=2.5,height=2.5)
  set_margin()

  m = 0;
  for(j in 1:nfiles) { m = max(m, as.numeric(data[[j]][field,])); }

  plot(kmers,as.numeric(data[[1]][field,]),type='b',
       col=cols[1], pch=pchtypes[1], lty=linetypes[1],
       ylim=c(0,2*m),axes=F,ylab=description,xlab="kmer",
       )#main=paste(description,"vs kmer-size"))
  axis(side = 2)
  axis(side = 1,at=kmers)

  for(j in 2:nfiles) {
    points(kmers,as.numeric(data[[j]][field,]),type='b',
           col=cols[j], pch=pchtypes[j], lty=linetypes[j])
  }

  s0=0; sn=0;
  for(j in 1:nfiles) {
    s0 = s0 + as.numeric(data[[j]][field,1]);
    sn = sn + as.numeric(data[[j]][field,length(kmers)-1]);
  }

  lpos='topright';
  # if(s0 < sn) { lpos='bottomright' }

  # legtxt = c('Perfect Plain', 'Perfect Links', 'Perfect String',
  #            'Stoch. Plain', 'Stoch. Links', 'Stoch. String',
  #            'Stoch. Err Plain', 'Stoch. Err Links', 'Stoch. Err String');
  legtxt=titles

  legend(lpos,legtxt,col=cols,pch=pchtypes,lty=linetypes)
}

dev.off()
