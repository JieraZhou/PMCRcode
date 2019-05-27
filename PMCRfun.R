### PMCR: function for Precision Medicine with Competing Risks (2 Med options, 2 risks)
library(rgenoud)
library(rpart)

# generate data
gen.data<-function(n,q=0.8,eta=c(0,1,-1,0.5),c0=8,truemod=TRUE){
  Z1<-runif(n,-2,2)
  Z2<-runif(n,-2,2)
  if(truemod) eze<-exp(cbind(1,Z1,Z2)%*%eta[-length(eta)])
  if(!truemod) eze<-exp(cbind(1,Z1,Z2,Z1^2)%*%eta)
  A<-rbinom(n,1,eze/(1+eze))
  
  eza1<-exp(-Z1+A*(Z1-Z2))      # exponential of linear predictors in PH model of risk 1
  eza2<-exp(-A*(2+Z1+2*Z2))     # rate of expoential dist. for risk2

  p1<-1-(1-q)^eza1
  rsk1<-rbinom(n,1,p1)

  U<-runif(n)
  T1<--log(1-(1-(1-p1*U)^(1/eza1))/q)
  T2<--log(1-U)/eza2
  T<-rsk1*T1+(1-rsk1)*T2
  ebs<-rsk1+2*as.numeric(rsk1==0)

  C<-runif(n,0,c0)
  obst<-apply(cbind(T,C),1,min)
  delta<-as.numeric(T<=C)
  status<-delta*ebs
return(data.frame(id=1:n,time=obst,status=status,Z1=Z1,Z2=Z2,A=A,piz=eze/(1+eze)))
}

# function calculating the CIFs
calc.F<-function(beta,mydata,pix.hat,Xp){
  xb<-Xp%*%beta
  n<-nrow(mydata)
  # smoothed version of weight
  h<-sd(xb)*(n/4)^(-1/3)
  w.numer<-pnorm(xb/h)*mydata$A+(1-pnorm(xb/h))*(1-mydata$A)
  w.denom<-pix.hat*mydata$A+(1-pix.hat)*(1-mydata$A)
  wi<-w.numer/w.denom

  ord.t0<-sort(mydata$time)
  ord.d0<-as.numeric(mydata$status>0)[order(mydata$time)]
  ord.ebs<-mydata$status[order(mydata$time)]
  ord.wi<-wi[order(mydata$time)]

  tt<-ord.t0[ord.d0==1]
  SIt<-numeric()
  SIt[1]<-1
  for(k in 1:length(tt)){
    SIt[k+1]<-SIt[k]*(1-sum(ord.wi[ord.t0==tt[k]])/sum(ord.wi[ord.t0>=tt[k]]))
  }

  tt1<-ord.t0[ord.ebs==1]
  tt2<-ord.t0[ord.ebs==2]
  lamb1<-lamb2<-numeric()
  for(k in 1:length(tt1)) lamb1[k]<-sum(ord.wi[ord.t0==tt1[k]])/sum(ord.wi[ord.t0>=tt1[k]])
  for(k in 1:length(tt2)) lamb2[k]<-sum(ord.wi[ord.t0==tt2[k]])/sum(ord.wi[ord.t0>=tt2[k]])

  if(length(tt1)>0){ F1.tt0<-stepfun(tt1,cumsum(c(0,stepfun(tt,SIt)(tt1)*lamb1)))}else{F1.tt0<-0}
  if(length(tt2)>0){ F2.tt0<-stepfun(tt2,cumsum(c(0,stepfun(tt,SIt)(tt2)*lamb2)))}else{F2.tt0<-0}
return(list(F1t0=F1.tt0,F2t0=F2.tt0,St0=stepfun(tt,SIt)))
}

# optim function
opt.beta1<-function(beta,mydata,pix.hat,Xp,tt0=5){
   cif.fit<-calc.F(beta,mydata,pix.hat,Xp)
   if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
return(F1.tt0)
}
opt.beta2<-function(beta,mydata,pix.hat,Xp,tt0=5){
   cif.fit<-calc.F(beta,mydata,pix.hat,Xp)
   if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
return(F2.tt0)
}
opt.beta<-function(beta,mydata,pix.hat,Xp,tt0=5,alp=0.1,M=1000){
   cif.fit<-calc.F(beta,mydata,pix.hat,Xp)
   if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
   if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
return(F1.tt0+M*(F2.tt0-alp)*as.numeric(F2.tt0>alp))
}

####################################################################################################
# Arguments in function PMCR                                                                       #
# Time: observed event or censoring time                                                           #
# Event: indicator variable of events or censoring(0)                                              #
# formula: treatment group indicator ~ predictors                                                  #
# data: data frame with all variables needed                                                       #
# rgenoud: if TRUE, the rgenoud function will be applied to search beta                            #
# Restrict: if TRUE, restricted regime is estimated. otherwise, unrestricted regimes are returned. #
# propscore: method for calculating the propensity scores. "logistic" or "tree".                   #
#            Note, the option "true" can be specified when the data is generated using gen.data,   #
#            and the true propensity score will be used.                                           #
# t0: time for CIF to be evaluated at.                                                             #
# alp: alpha value as the upper limit of risk2 CIF                                                 #
# M: a large number used as penalty for the restriction being violated                             #
####################################################################################################

PMCR<-function(Time,Event,formula=formula(data),data=parent.frame(),rgenoud=TRUE,Restrict=TRUE,propscore="logistic",t0=5,alp=0.1,M=1000){
  mydata<-data
  call <- match.call()
  avars<-all.vars(formula)
  mydata$time<-data[,Time]
  mydata$status<-data[,Event]
  mydata$A<-data[,avars[1]]
  if(propscore=="true")     mydata$pix.hat<-mydata$piz
  if(propscore=="logistic") mydata$pix.hat<-predict(glm(formula,mydata,family=binomial("logit")),type="response")
  if(propscore=="tree")     mydata$pix.hat<-predict(rpart(formula,mydata,method="class"))[,"1"]
  temp.x<-terms(formula, data = mydata)
  Xp <- model.matrix(temp.x,mydata)
  colnames(Xp)<-paste("Z",c(1:ncol(Xp))-1,sep="")

  npar<-ncol(Xp)
  beta0<-matrix(0,2*npar,npar)
  for(i in 1:npar){
    beta0[2*i-1,i]<-1
    beta0[2*i,i]<--1
  }
  outest<-list(data=mydata)

  if(!Restrict){
    f1val<-1
    for(k in 1:nrow(beta0)){
      fit1<-try(optim(par=beta0[k,],fn=opt.beta1,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
      if(fit1$value<f1val){
        outest$beta1<-fit1$par/sqrt(sum(fit1$par^2))
        f1val<-fit1$value
        #cat(paste("k=",k,"\n"))
        #print(c(outest$beta1,f1val))
      }
    }

    f2val<-1
    for(k in 1:nrow(beta0)){
      fit2<-try(optim(par=beta0[k,],fn=opt.beta2,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
      if(fit2$value<f2val){
        outest$beta2<-fit2$par/sqrt(sum(fit2$par^2))
        f2val<-fit2$value
        #cat(paste("k=",k,"\n"))
        #print(c(outest$beta2,f2val))
      }
    }
    

  if(rgenoud){
    fit1<-try(genoud(opt.beta1,nvars=npar,max=FALSE,starting.values=fit1$par,max.generations=30,print.level=0,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
    if(!is.character(fit1)){
      if(fit1$value<f1val){
        outest$beta1<-fit1$par/sqrt(sum(fit1$par^2))
        f1val<-fit1$value
        #cat("rgenoud replace\n")
        #print(c(outest$beta1,f1val))
      }
    }
    fit2<-try(genoud(opt.beta2,nvars=npar,max=FALSE,starting.values=fit2$par,max.generations=30,print.level=0,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
    if(!is.character(fit2)){
      if(fit2$value<f2val){
        outest$beta2<-fit2$par/sqrt(sum(fit2$par^2))
        f2val<-fit2$value
        #cat("rgenoud replace\n")
        #print(c(outest$beta2,f2val))
      }
    }
  }
  FF<-calc.F(beta=outest$beta1,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp)
  outest$Fbeta1<-c(FF$F1t0(t0),FF$F2t0(t0))
  outest$f1val<-f1val
  FF<-calc.F(beta=outest$beta2,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp)
  outest$Fbeta2<-c(FF$F1t0(t0),FF$F2t0(t0))
  outest$f2val<-f2val
  }

  if(Restrict){
    f3val<-1
    for(k in 1:nrow(beta0)){
      fit<-try(optim(par=beta0[k,],fn=opt.beta,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0,alp=alp,M=M),silent=TRUE)
      if(fit$value<f3val){
        outest$beta3<-fit$par/sqrt(sum(fit$par^2))
        f3val<-fit$value
        #cat(paste("k=",k,"\n"))
        #print(c(outest$beta3,f3val))
      }
    }
    if(rgenoud){
      fit<-try(genoud(opt.beta,nvars=npar,max=FALSE,starting.values=fit$par,max.generations=30,print.level=0,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0,alp=alp,M=M),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<f3val){
          outest$beta3<-fit$par/sqrt(sum(fit$par^2))
          f3val<-fit$value
          #cat("rgenoud replace\n")
          #print(c(outest$beta3,f3val))
        }
      }
    }
    FF<-calc.F(beta=outest$beta3,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp)
    outest$f3val<-f3val
    outest$Fbeta3<-c(FF$F1t0(t0),FF$F2t0(t0))
  }
return(outest)
}