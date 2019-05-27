
### fit for the pseudoHIV data ###
source("PMCRfun.R")

pseudoHIV<-read.table("pseudoHIV.txt",head=TRUE)
head(pseudoHIV)

t0<-365
# Unrestricted regime
temp.fit1<-PMCR(Time="time",Event="status",formula=A~Age+Gender+Black+GC,data=pseudoHIV,rgenoud=FALSE,Restrict=FALSE,propscore="logistic",t0=t0)
temp.fit1$beta1  #unrestricted regime for risk1
temp.fit1$Fbeta1 #CIFs for risk1 and risk2 with regime beta1
temp.fit1$beta2  #unrestricted regime for risk1
temp.fit1$Fbeta2 #CIFs for risk1 and risk2 with regime beta2

#range of alpha
alps<-c(temp.fit1$Fbeta2[2],seq(round(temp.fit1$Fbeta2[2],2)+0.01,round(temp.fit1$Fbeta1[2],2)+0.03,0.01))

#Restricted regime
alp<-0.1;M<-1e+5
temp.fit2<-PMCR(Time="time",Event="status",formula=A~Age+Gender+Black+GC,data=pseudoHIV,rgenoud=FALSE,Restrict=TRUE,propscore="logistic",t0=t0,alp=alp,M=M)
temp.fit2$beta3
temp.fit2$Fbeta3
