
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
  stop("Usage: R -f make_link_cutoffs.R --args <kmer> <links.csv>");
}

kmer<-as.integer(args[1])
links_csv<-args[2]

fit_gamma_poiss_model <- function(data,stepsize=0.005)
{
  data<-data+0.0001 # To avoid div by zero
  r1<-data[2]/data[1];
  r2<-data[3]/data[2];
  rr<-r2/r1;
  aa<-seq(stepsize,2,stepsize);
  f.aa<-gamma(aa)*gamma(aa+2)/(2*gamma(aa+1)^2);
  a.est<-aa[which.min(abs(f.aa-rr))];
  b.est<-gamma(a.est+1)/(r1*gamma(a.est))-1;
  b.est<-max(b.est,0.000001)

  c0<-data[1]*(b.est/(1+b.est))^(-a.est);

  aa<-1:length(data);
  e.cov<-a.est*log(b.est)-lgamma(a.est)-lfactorial(aa-1)+lgamma(a.est+aa-1)-(a.est+aa-1)*log(1+b.est);
  e.cov<-exp(e.cov)*c0;

  return(e.cov)
}

# Links shorter than kmer+2 shouldn't exist mmkay
maxcount <- 150
minlen <- kmer+2
maxlen <- kmer+10
krange <- minlen:maxlen;
thresh <- 0.05

d <- read.table(links_csv, as.is=T, head=T, sep=",");

if(nrow(d) == 0) {
  warning(paste("No links in",links_csv));
  cat(0,'\n');
  quit(save="no",status=0);
}

z <- matrix(0,maxlen,maxcount)
cutoffs <- c()

# Construct histograms
for(linklen in krange) {
  c <- d[d[,'LinkLength']==linklen,]
  h <- hist(c[c[,'Count']<=maxcount,'Count'],breaks=0:maxcount,plot=F)
  z[linklen,h$mids+0.5] = h$counts
  z[linklen,maxcount] = z[linklen,maxcount] + sum(c[c[,'Count']>maxcount,'Count'])
}

for(linklen in krange) {
  e.cov <- fit_gamma_poiss_model(z[linklen,])
  data <- pmax(z[linklen,]-e.cov,1)
  cut <- 0
  for(i in 3:maxcount) {
    if(e.cov[i]/(data[i]+e.cov[i]) < thresh) {
      cut <- i
      break
    }
  }
  if(cut == 0) { warn("Couldn't find cut"); }
  else { cutoffs[length(cutoffs)+1] <- cut; }
}

cat(paste(cutoffs),'\n');
cat(median(cutoffs),'\n');
