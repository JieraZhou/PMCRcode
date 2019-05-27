########################################
# Change the following directory w.r.t. the code and results folder
source("read_fun.R")

n=500                          #sample size
# table 2
q<-0.2;alps<-c(0.3,0.4,0.5)
#c0=4 ;RC=40   #40% rc, 28% risk 1, 32% risk 2
c0=55;RC=15   #15% rc, 35% risk 1, 50% risk 2
truth<-matrix(c(0.661,0.302,0,-0.05,-0.53,-0.707,0.748,0.792,0.707,0.163,0.137,0.134),3,4)
# table 3
#q<-0.5;alps<-c(0.2,0.3,0.4)
#c0=2.5;RC=40  #40% rc, 47% risk 1, 13% risk 2
#c0=28 ;RC=15  #15% rc, 62% risk 1, 23% risk 2
#truth<-matrix(c(0.609,0.163,1,-0.096,-0.604,-0.707,0.787,0.78,0.707,0.353,0.31,0.307),3,4)
# table 4
#q<-0.8;alps<-c(0.05,0.15,0.25)
#c0=1.8;RC=40  #40% rc, 55% risk 1,  5% risk 2
#c0=12 ;RC=15  #15% rc, 76% risk 1,  9% risk 2
#truth<-matrix(c(0.728,0.374,0,0.272,-0.336,-0.707,0.629,0.864,0.707,0.613,0.497,0.475),3,4)
################
eta0<-c(0,1,-0.5,0.5)               #coefficients in the logistic regression for propensity score

set.seed(100)
nsimu=1000
mydata<-list()
for(simu in 1:nsimu){
mydata[[simu]]<-gen.data(n=n,q=q,eta=eta0,c0=c0)
}
##############
#change the prefix names to match saved results in different settings
prefix.true<-paste("Case1_true_n",n,"ces",RC,"q",q,sep="") 
col.true<-resltabl(q=q,alps=alps,prefix=prefix.true,mydata=mydata,truth=truth)
prefix.logit<-paste("Case1_logit_n",n,"ces",RC,"q",q,sep="") 
col.logit<-resltabl(q=q,alps=alps,prefix=prefix.logit,mydata=mydata,truth=truth)
prefix.tree<-paste("Case1_tree_n",n,"ces",RC,"q",q,sep="") 
col.tree<-resltabl(q=q,alps=alps,prefix=prefix.tree,mydata=mydata,truth=truth)


