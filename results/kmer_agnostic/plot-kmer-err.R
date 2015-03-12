# R script for plotting experimental results

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

set_margin <- function() {
  par(oma=c(0,0,0,0))
  par(mar=c(5,4,1.5,1))
}

args <- commandArgs(trailingOnly = TRUE)

# files=c()
# kmers= c(15,21,31,41,51,63,75,99)
# for(i in 1:length(kmers)) { files[i] = paste('k',kmers[i],'/kmer_cleaning/kmer_cleaning.csv',sep='') }

out_path=args[1]
files=args[2:length(args)]
nfiles=length(files)
titles=substr(files,1,3)

# Thresholds used for cleaning, pulled from k*/graphs/stocherr.ctx.log:
# for k in k*; do cat $k/graphs/stocherr.ctx.log | grep 'cleaning threshold'; done | grep -o '[0-9]*$'
threshs=c(9,8,8,9,8,7,5)

cols=rainbow(nfiles);

data=list()
maxErrKmers=0
for(i in 1:nfiles) {
  data[[i]]=read.csv(file=files[i],as.is=T,header=T)
  cat("files[",i,"] = '",files[i],"'\n",sep='')
  maxErrKmers=max(maxErrKmers,max(data[[i]][,'errorKmers']))

  t=threshs[i]
  x=which.min(abs(data[[i]][,'kmerThresh']-t))
  act=data[[i]][x,'kmerThresh']
  cat(titles[i],' thresh=',t,'; at ',act,' kmer mismatch = ',data[[i]][x,'errorKmers']/(data[[i]][x,'correctKmers']+data[[i]][x,'errorKmers']),'\n')
}

pdf(file=out_path,width=10,height=5)

par(mfrow=c(1,2))

i=1
plot(data[[i]][,'kmerThresh'], 100*data[[i]][,'errorKmers']/data[[i]][,'correctKmers'],
     xlim=c(0,30), ylim=c(0,10),
     col=cols[i],
     xlab='Kmer Threshold', ylab='Match Rate (%)',
     main="Kmer kmerThresh vs matching rate",type='b');

for(i in 2:nfiles) {
  points(data[[i]][,'kmerThresh'], data[[i]][,'errorKmers']/data[[i]][,'correctKmers'], type='b', col=cols[i]);
}

legend('bottomright',titles,fill=cols)

i=1
plot(data[[i]][,'kmerThresh'], data[[i]][,'errorKmers'],
     xlim=c(0,30), ylim=c(0,maxErrKmers),
     xlab='Kmer Threshold', ylab='Number of bad kmers',
     main="kmer-thresh vs # bad kmers",type='b');

for(i in 2:nfiles) {
  points(data[[i]][,'kmerThresh'],data[[i]][,'errorKmers'], type='b', col=cols[i]);
}

legend('topright',titles,fill=cols)

dev.off();

