## install.packages
if (!require('data.table')) install.packages('data.table')
if (!require('dplyr')) install.packages('dplyr')
if (!require('BGLR')) install.packages('BGLR')

## phenotype
ref_precor_phe = fread("ref_precor_phe.txt",data.table=F)

## genotype
geno=fread("mice_835_qc_GRCm39.raw",data.table=F)

n=nrow(geno)
p=ncol(geno)-6
X=geno[,7:ncol(geno)]

## Centering and standardization
for(i in 1:p)
{ 
  X[,i]<-(X[,i]-mean(X[,i]))/sd(X[,i]) 
}

## bayesian method
y=ref_precor_phe$bmd1

nIter=120000;
burnIn=20000;
thin=100;
saveAt='';
S0=NULL;
weights=NULL;
R2=0.5;
ETA<-list(list(X=X,model='BayesB',probIn=0.05))

fit=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)

save(fit, file = "bb_bmd1/fit_BB.RData")
