
# R --vanilla -f plot-results.R --args <links.csv> <plain.csv>
#options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  message("Usage: R --vanilla -f plot-confidence.R --args <out.pdf>")
  quit()
}

outpath=args[1]

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

contig_conf <- function(covg,readlen,kmer_size)
{
  if(kmer_size > readlen) { return(0); }
  lambda = covg / readlen;
  read_kmers = readlen - kmer_size + 1;
  power = (1 - exp(-lambda * read_kmers)) *
          exp(-lambda * exp(-lambda * read_kmers));
  return(power);
}

covgs_a=c(10,10,30,60,100)
rlens_a=c(75,100,100,100,100)

covgs_b=c(30,30,60,100)
rlens_b=c(150,250,250,250)

# covgs=c(covgs_a,covgs_b)
# rlens=c(rlens_a,rlens_b)

lens=c(75,100,150,250)
rlens=rep(lens,each=10)
covgs=rep(1:10*10,4)

N=length(covgs)

cols=rainbow(N)

m=max(rlens)+10
conf=c()
txt=c()
linetypes=c()

pdf(file=outpath,width=6,height=6)

plot(NA,xlim=c(0,m),ylim=c(0,1),xlab="Distance between junctions (bp)",ylab="Confidence",
     main="Confidence in assembling between junctions\nfor given read lengths and coverage 10,20,...,100X")

for(i in 1:N) {
  covg=covgs[i]
  rlen=rlens[i]
  for(j in 1:(rlen+1)) { conf[j] = contig_conf(covg,rlen,j) }
  # linetypes[i]=((i-1) %% 6) + 1
  linetypes[i]=1
  # col=cols[floor((i-1)/10)+1]
  col=cols[i]
  points(1:(rlen+1), conf[1:(rlen+1)], type='l',lty=linetypes[i],col=col);
  txt[i] = paste(covg,"X ",rlen,"bp",sep='')
}

# legend('bottomleft',txt,lty=linetypes)
legend('bottomleft',paste(lens,'bp',sep=''),fill=cols[c(5,15,25,35)],
       title="Read lengths")

dev.off()
