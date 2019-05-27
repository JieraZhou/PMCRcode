source("PMCRfun.R")

# Simulation Parameters
n=500                    #sample size
eta0<-c(0,1,-0.5,1)      #coefficients in the logistic regression for propensity score

set.seed(100)
nsimu=3                  #number of replications        

# Similation Settings
#alps: list of alpha values used for constrain the risk2 CIF.
#c0: constant used to adjust the right censoring proportions.

# simulation table 2
#q<-0.2;alps<-c(0.3,0.4,0.5)
#c0=4 ;RC=40   #40% rc, 28% risk 1, 32% risk 2
#c0=55;RC=15   #15% rc, 35% risk 1, 50% risk 2

# simulation table 3
#q<-0.5;alps<-c(0.2,0.3,0.4)
#c0=2.5;RC=40  #40% rc, 47% risk 1, 13% risk 2
#c0=28 ;RC=15  #15% rc, 62% risk 1, 23% risk 2

# simulation table 4
q<-0.8;alps<-c(0.05,0.15,0.25)
#c0=1.8;RC=40  #40% rc, 55% risk 1,  5% risk 2
c0=12 ;RC=15  #15% rc, 76% risk 1,  9% risk 2

###################################################################################################################
# Instructions:                                                                                                   #
# the following specification corresponds to the Case 1 & logistic propensity score in simulation                 #
# To derive other results:                                                                                        #
# (1) adjust "truemod" in gen.data to have different models: truemod=TRUE for Case 1 and truemod=FALSE for Case 2 #
# (2) change "propscore" in fit1-fit3 to specify different propensity scores: "true", "logistic" or "tree"        #
###################################################################################################################
#setwd("C:/Users/gazhou/Dropbox/Research/Personalized Medicine/paper/JASA/submit/Code/results")
prefix<-paste("Case1_logit_n",n,"ces",RC,"q",q,sep="") #change the prefix names when saving results for different settings
tt0<-2;M<-1e+5
beta.est1<-beta.est2<-beta.est3<-matrix(nrow=nsimu,ncol=3)
F1.est<-F2.est<-fval<-attained<-matrix(nrow=nsimu,ncol=3)

for(simu in 1:nsimu){
  mydata<-gen.data(n=n,q=q,eta=eta0,c0=c0,truemod=TRUE)  

  ptm<-proc.time()
    fit1<-try(PMCR(Time="time",Event="status",A~Z1+Z2,data=mydata,rgenoud=FALSE,Restrict=TRUE,propscore="logistic",t0=tt0,alp=alps[1],M=M),silent=TRUE)
  rectm1<-proc.time()-ptm

  ptm<-proc.time()
    fit2<-try(PMCR(Time="time",Event="status",A~Z1+Z2,data=mydata,rgenoud=FALSE,Restrict=TRUE,propscore="logistic",t0=tt0,alp=alps[2],M=M),silent=TRUE)
   rectm2<-proc.time()-ptm

  ptm<-proc.time()
    fit3<-try(PMCR(Time="time",Event="status",A~Z1+Z2,data=mydata,rgenoud=FALSE,Restrict=TRUE,propscore="logistic",t0=tt0,alp=alps[3],M=M),silent=TRUE)
  rectm3<-proc.time()-ptm

  if(!is.character(fit1) & !is.character(fit2) &!is.character(fit3)){
    beta.est1[simu,]<-fit1$beta3
    beta.est2[simu,]<-fit2$beta3
    beta.est3[simu,]<-fit3$beta3

    fval[simu,]<-c(fit1$f3val,fit2$f3val,fit3$f3val)
    F1.est[simu,]<-c(fit1$Fbeta3[1],fit2$Fbeta3[1],fit3$Fbeta3[1])
    F2.est[simu,]<-c(fit1$Fbeta3[2],fit2$Fbeta3[2],fit3$Fbeta3[2])
    attained[simu,]<-as.numeric(F2.est[simu,]<=alps)

    cens<-c(mean(mydata$status==0),mean(mydata$status==1),mean(mydata$status==2))

    # output and save the estimated results item-by-item
    cat("simu=",simu,",cen=",cens[1],"\n")
    print(cbind(beta.est1[simu,],beta.est2[simu,],beta.est3[simu,]))
    cat("attained=",F2.est[simu]<=alps,"\n")
    write.table(t(c(beta.est1[simu,],beta.est2[simu,],beta.est3[simu,])),file=paste(prefix,"_beta.txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(fval[simu,]),file=paste(prefix,"_fval.txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(F1.est[simu,]),file=paste(prefix,"_F1est.txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(F2.est[simu,]),file=paste(prefix,"_F2est.txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(attained[simu,]),file=paste(prefix,"_attain.txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(matrix(c(rectm1[1],rectm2[1],rectm3[1]),1,3),paste(prefix,"rectm.txt",sep=""),row.names=FALSE,col.names=FALSE,append=TRUE)
    write.table(t(cens),file=paste(prefix,"_prop.txt",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
  }
}