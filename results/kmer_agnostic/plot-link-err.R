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
# kmers= c(15,21,31,41,51,63,75)
# for(i in 1:length(kmers)) { files[i] = paste('k',kmers[i],'/link_cleaning/link_cleaning.csv',sep='') }

out_path=args[1]
files=args[2:length(args)]
nfiles=length(files)
titles=substr(files,1,3)

# Thresholds used for cleaning, pulled from:
# for k in k*; do echo $k `tail -1 $k/links/cleaning.txt`; done
threshs=c(17,16,6,5,4,3,3,0)

cols=rainbow(nfiles);

data=list()
maxNumLinks=0
for(i in 1:nfiles) {
  data[[i]]=read.csv(file=files[i],as.is=T,header=T)
  cat("files[",i,"] = '",files[i],"'\n",sep='')
  maxNumLinks=max(maxNumLinks,max(data[[i]][,'numLinks']))

  t=threshs[i]
  x=which.min(abs(data[[i]][,'linkThresh']-t))
  act=data[[i]][x,'linkThresh']
  cat(titles[i],' thresh=',t,'; at ',act,' link mismatch = ',1-data[[i]][x,'numMatch']/data[[i]][x,'numLinks'],'\n')
}

pdf(file=out_path,width=10,height=5)

par(mfrow=c(1,2))

i=1
plot(data[[i]][,'linkThresh'],data[[i]][,'matchRate'],ylim=c(0,100),
     col=cols[i],
     xlab='Link Threshold', ylab='Match Rate (%)',
     main="Link threshold vs matching rate",type='b');

for(i in 2:nfiles) {
  points(data[[i]][,'linkThresh'],data[[i]][,'matchRate'], type='b', col=cols[i]);
}

legend('bottomright',titles,fill=cols)

i=1
plot(data[[i]][,'linkThresh'],data[[i]][,'numLinks'],ylim=c(0,maxNumLinks),
     xlab='Link Threshold', ylab='Number of links',
     main="Link threshold vs number of links",type='b');

for(i in 2:nfiles) {
  points(data[[i]][,'linkThresh'],data[[i]][,'numLinks'], type='b', col=cols[i]);
}

legend('topright',titles,fill=cols)

dev.off();

