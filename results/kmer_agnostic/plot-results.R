# R script for plotting experimental results

# R --vanilla -f plot-results.R --args <links.csv> <plain.csv>
options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  message("Usage: R --vanilla -f plot-results.R --args <links.csv> <plain.csv>")
  quit()
}

links_file=args[1]
plain_file=args[2]

links=read.csv(file=links_file,as.is=T,row.names=1)
plain=read.csv(file=plain_file,as.is=T,row.names=1)

kmers=as.numeric(links['kmer',])
cols=c('firebrick','dodgerblue')

fields=c('median','N50','length','contigs')
descriptions=c("Median contig length","Contig N50","Assembled length","Number of contigs")
path='plots/'

for(i in 1:length(fields)) {
  field = fields[i]
  description = descriptions[i]
  pdf(file=paste(path,field,'.pdf',sep=''),width=7,height=7)
  m=max(as.numeric(links[field,]),as.numeric(plain[field,]))
  plot(kmers,as.numeric(links[field,]),type='b',col=cols[1],
       ylim=c(0,m),axes=F,ylab=description,xlab="kmer",
       main=paste(description,"vs kmer-size"))
  axis(side = 2)
  axis(side = 1,at=kmers)
  points(kmers,as.numeric(plain[field,]),type='b',col=cols[2])
  lpos='topright';
  if(plain[field,length(kmers)-1] > plain[field,1]) { lpos='bottomright' }
  legend(lpos,c('Links','Plain'),fill=cols)
  dev.off()

  # Plot links/plain separately
  for(j in 1:2) {
    if(j == 1) { src_str = 'links'; src=links; }
    else       { src_str = 'plain'; src=plain; }
    m=max(as.numeric(src[field,]))
    pdf(file=paste(path,field,'.',src_str,'.pdf',sep=''),width=7,height=7)
    plot(kmers,as.numeric(src[field,]),type='b',col=cols[j],
         ylim=c(0,m),axes=F,ylab=description,xlab="kmer",
         main=paste(description," vs kmer-size (",src_str," only)",sep=''))
    axis(side = 2)
    axis(side = 1,at=kmers)
    dev.off()
  }
}
