

#To calculate error threshold

#data<-read.table("ERR019054.unclean.kmer31.q10.ctx.covg", as.is=T, head=T);
data<-read.table("na12878_exome.k31.q10.ctx.covg", as.is=T, head=T);

do.sim<-FALSE;
if (do.sim) {
#Simulate data
a.sim<-0.5;
b.sim<-1/5;
vv<-rgamma(1e5, a.sim, rate=b.sim);
dd<-rpois(1e6, vv);
tb<-table(dd);
data<-array(0, c(100,2));
data[,1]<-0:(nrow(data)-1);
data[,2]<-c(0, tb[1:99]);
}

##########################################################################################
##Version 1 - assumes all kmers with coverage 1 or 2 are error and fits exponential decay
##########################################################################################

#Fit exponential decay model
ll<-log(data[2,2]/data[3,2]);

#Plot
plot(data, log="xy");

n.max<-100; 	#Assumes everything above n.max is signal
aa<-1:n.max;
e.cov<-data[2,2]*exp((aa-1)*(-ll));
points(aa, e.cov, col="red", pch=19);

fdr<-1-(data[(1:n.max)+1,2]-e.cov)/data[(1:n.max)+1,2];
thresh<-0.999;

pt<-which(fdr>(1-thresh));
cut<-max(pt)+1;

abline(v=cut, col="blue");

#For P. falciparum data indicates that coverage <=5 should be treated as error


##########################################################################################
#Version 2 - assumes all kmers with coverage <=3 are error and fits gamma model
##########################################################################################

#Fit gamma decay model
r1<-data[3,2]/data[2,2];
r2<-data[4,2]/data[3,2];
rr<-r2/r1;

aa<-seq(0.01,2,0.01);
f.aa<-gamma(aa)*gamma(aa+2)/(2*gamma(aa+1)^2); # r2/r1
# Estimate alpha, beta for gamma distribution
a.est<-aa[which.min(abs(f.aa-rr))];
b.est<-gamma(a.est+1)/(r1*gamma(a.est))-1;
c0<-data[2,2]*(b.est/(1+b.est))^(-a.est);

#Plot
plot(data, log="xy");

n.max<-100; 	#Assumes everything above n.max is signal
aa<-1:n.max;
e.cov<-a.est*log(b.est)-lgamma(a.est)-lfactorial(aa-1)+lgamma(a.est+aa-1)-(a.est+aa-1)*log(1+b.est);
e.cov<-exp(e.cov)*c0;
points(aa, e.cov, col="red", pch=19);

fdr<-1-(data[(1:n.max)+1,2]-e.cov)/data[(1:n.max)+1,2];
thresh<-0.999;

pt<-which(fdr>(1-thresh));
cut<-max(pt)+1;

abline(v=cut, col="blue");

#For P. falciparum data indicates that coverage <=10 should be treated as error


