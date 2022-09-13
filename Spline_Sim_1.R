#####################################################################
#                      Load Packages/Data Generation
#####################################################################

library(ECIC1.1Dev)
library(paleoTS)
library(splines)

#####################################################################
#                      Define Functions
#####################################################################
# This function computes the KL divergence for two MVN models
KLforMVN <- function(Sig1,Sig2,n,mu1,mu2) #Cov Mat's, dimension of MVN model, and means 
{
  invSig2 <- solve(Sig2)
  logDets<- log(det(Sig2)/det(Sig1))
  sigTr <- sum(diag(invSig2%*%Sig1))
  musPart <- t(mu2-mu1)%*%invSig2%*%(mu2-mu1)
  KLDiv <- 0.5*(logDets-n+sigTr+musPart)
  return(KLDiv)
}

# intermediary function for pseudo true coefficients
SigHat <- function(x,w=NULL)
{
  #store the number of observations
  n <- nrow(x)
  if(is.null(w))
  {
    SigEst <- x%*%t(x) 
  }
  else
  {
    SigEst <- x%*%t(w)
  }
  return(SigEst/n)
}
# compute the psuedo-true parameter values as in pg. 7 in 
# Pesaran & Weeks
psuedoTrueCoefs <- function(trueCoef,trueVar,trueBases,altBases)
{
  SigxxInv <- solve(SigHat(altBases))
  Sigww <- SigHat(trueBases)
  Sigwx <- SigHat(trueBases,altBases)
  Sigxw <- SigHat(altBases,trueBases)
  coefEst <- SigxxInv%*%Sigxw%*%trueCoef
  varEst <- trueVar + t(trueCoef)%*%(
    Sigww-Sigwx%*%SigxxInv%*%Sigxw
  )%*%trueCoef
  paramEst <- list(coefs=coefEst,var=as.numeric(varEst))
  return(paramEst)
}
#####################################################################
#                      Generate different spline fits
#####################################################################

# where to store spline fits
# wrong spline fit no1 (W for wrong)
splinesW1 <- list()
dfW1 <- 3
# wrong spline fit no2 (W for wrong)
splinesW2 <- list()
dfW2 <- 5
# wrong spline fit, but closest to correct  (Cl for close)
splinesCl <- list()
dfCl <- 8
# spline data that correctly smooths (Co for correct)
splinesCo <- list()
# set the true coefficients and variance for the spline
# knots are at endpoints and quintiles
coefCo <- c(1,-1,1,-1,1,-1,1)
varCo <- 0.0036
dfCo <- 7
# list to store observations
ys <- list()

# store the pseudo true parameters for the wrongs fits
pseudParsW1 <- list()
pseudParsW2 <- list()
pseudParsCl <- list()

# randomly sample values between -50 and 50 for the x values
lowLim <- -10
upLim <- 10
# set three different total observation points
ns <- c(20) #50, 100,400,800
# store the observation points in a list
xs <- list()
KLMat <- matrix(NA,nrow=3,ncol=length(ns))
colnames(KLMat) <- ns
rownames(KLMat) <- c("W1","W2","Cl")
manyX <- seq(from=lowLim,to=upLim,length.out=200)
manyBase <- bs(manyX,df=dfCo,degree=3)
manyY <- manyBase%*%coefCo
for(i in 1:length(ns))
{
  xs[[i]] <- seq(from=lowLim,to=upLim,length.out=ns[i])
  w1BaseTemp <- bs(xs[[i]],df=dfW1,degree=3)
  w2BaseTemp <- bs(xs[[i]],df=dfW2,degree=3)
  clBaseTemp <- bs(xs[[i]],df=dfCl,degree=3)
  coBaseTemp <- bs(xs[[i]],df=dfCo,degree=3)
  w1PseudParsTemp <- psuedoTrueCoefs(trueCoef=coefCo,trueVar=varCo,
                                    trueBases=t(coBaseTemp),altBases=t(w1BaseTemp))
  w2PseudParsTemp <- psuedoTrueCoefs(trueCoef=coefCo,trueVar=varCo,
                                     trueBases=t(coBaseTemp),altBases=t(w2BaseTemp))
  clPseudParsTemp <- psuedoTrueCoefs(trueCoef=coefCo,trueVar=varCo,
                                     trueBases=t(coBaseTemp),altBases=t(clBaseTemp))
  pseudParsW1[[i]] <- w1PseudParsTemp
  pseudParsW2[[i]] <- w2PseudParsTemp
  pseudParsCl[[i]] <- clPseudParsTemp
  ys[[i]] <- coBaseTemp%*%coefCo
  # compute the KL divergences for the wrong fits vs the true fit
  # using the quasi-ML estimators in Pesaran & Weeks 
  KLEstW1 <- KLforMVN(Sig1=diag(ns[i])*varCo,Sig2=diag(ns[i])*pseudParsW1[[i]][[2]],n=ns[i],
                     mu1=coBaseTemp%*%coefCo,mu2=w1BaseTemp%*%pseudParsW1[[i]][[1]])
  KLEstW2 <- KLforMVN(Sig1=diag(ns[i])*varCo,Sig2=diag(ns[i])*pseudParsW2[[i]][[2]],n=ns[i],
                      mu1=coBaseTemp%*%coefCo,mu2=w2BaseTemp%*%pseudParsW2[[i]][[1]])
  KLEstCl <- KLforMVN(Sig1=diag(ns[i])*varCo,Sig2=diag(ns[i])*pseudParsCl[[i]][[2]],n=ns[i],
                      mu1=coBaseTemp%*%coefCo,mu2=clBaseTemp%*%pseudParsCl[[i]][[1]])
  KLMat[,i] <- c(KLEstW1,KLEstW2,KLEstCl)
  # plot the pseudo-ML estimated means for the over and under smoothed model for each sample size
  plot(manyX,manyY,type="l",main=paste("Fits W/ Psuedo-ML Ests.,Sample Size=",ns[[i]],"\n
       Black=Co,Red=W1,Blue=W2,Purple=Cl"),
       xlab="x",ylab="y",ylim=c(min(ys[[i]])-0.2,max(ys[[i]])+0.2))
  lines(xs[[i]],w1BaseTemp%*%pseudParsW1[[i]][[1]],col="red")
  lines(xs[[i]],w2BaseTemp%*%pseudParsW2[[i]][[1]],col="blue")
  lines(xs[[i]],clBaseTemp%*%pseudParsCl[[i]][[1]],col="purple")
}

# generate 2,000 draws of each sample size 
set.seed(222)
drawSize <- 2000
datList <- list()
for(i in 1:length(ns))
{
  tempN <- ns[i]
  tempMat <- matrix(data=NA,nrow=tempN,ncol=drawSize)
  for(j in 1:drawSize)
  {
    # draw from a normal dist with mean 0 and sd .06
    tempDraw <- ys[[i]] + rnorm(n=tempN,mean=0,sd=0.06)
    tempMat[,j] <- tempDraw
  }
  datList[[i]] <- tempMat
}

# plot just the first draws for each sample size to illustrate what's going on
for(i in 1:length(ns))
{
  splinesW1[[i]] <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfW1,degree=3))
  splinesW2[[i]] <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfW2,degree=3))
  splinesCl[[i]] <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfCl,degree=3))
  splinesCo[[i]] <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfCo,degree=3))
  plot(xs[[length(xs)]],ys[[length(ys)]],type="l",main=paste("Spline Fits First Draw, Sample Size=",ns[[i]],"\n
       Points=Data,Black=Co,Orange=Co Fit,Red=W1,Blue=W2,Purple=Cl"),
      xlab="x",ylab="y",ylim=c(min(ys[[i]])-0.2,max(ys[[i]])+0.2))
  points(xs[[i]],datList[[i]][,1],pch=16,cex=0.5)
  lines(xs[[i]],splinesW1[[i]]$fitted.values,col="red")
  lines(xs[[i]],splinesW2[[i]]$fitted.values,col="blue")
  lines(xs[[i]],splinesCl[[i]]$fitted.values,col="purple")
  lines(xs[[i]],splinesCo[[i]]$fitted.values,col="darkorange3")
}

# determine which model is best by BIC and AIC with the following two model sets
# M1={W1,W2,Co} # Co should be chosen
# M2={W1,W2,Cl} # Cl should be chosen
# compute BICs and AICs for each simulated data sets
AICList <- list()
BICList <- list()
for(i in 1:length(ns))
{
  splinesTempW1 <- apply(X=datList[[i]],MARGIN=2,FUN=function(x) lm(x~bs(xs[[i]],df=dfW1,degree=3)))
  splinesTempW2 <- apply(X=datList[[i]],MARGIN=2,FUN=function(x) lm(x~bs(xs[[i]],df=dfW2,degree=3)))
  splinesTempCl <- apply(X=datList[[i]],MARGIN=2,FUN=function(x) lm(x~bs(xs[[i]],df=dfCl,degree=3)))
  splinesTempCo <- apply(X=datList[[i]],MARGIN=2,FUN=function(x) lm(x~bs(xs[[i]],df=dfCo,degree=3)))
  w1BICs <- sapply(X=splinesTempW1,FUN=function(x) BIC(x))
  w2BICs <- sapply(X=splinesTempW2,FUN=function(x) BIC(x))
  clBICs <- sapply(X=splinesTempCl,FUN=function(x) BIC(x))
  coBICs <- sapply(X=splinesTempCo,FUN=function(x) BIC(x))
  w1AICs <- sapply(X=splinesTempW1,FUN=function(x) AIC(x))
  w2AICs <- sapply(X=splinesTempW2,FUN=function(x) AIC(x))
  clAICs <- sapply(X=splinesTempCl,FUN=function(x) AIC(x))
  coAICs <- sapply(X=splinesTempCo,FUN=function(x) AIC(x))
  AICList[[i]] <- cbind(w1=w1AICs,w2=w2AICs,cl=clAICs,co=coAICs)
  BICList[[i]] <- cbind(w1=w1BICs,w2=w2BICs,cl=clBICs,co=coBICs)
  print(i)
}
names(AICList) <- names(BICList) <- paste("Sample Size=",ns,sep="")

# select the best model in both M1 and M2 strictly based on AIC and BIC
M1SelectAIC <- M2SelectAIC <- rep(NA,length(AICList))
M1SelectBIC <- M2SelectBIC <- rep(NA,length(BICList))
for(i in 1:length(BICList))
{
  # choose model based on BIC & AIC alone for # M1={W1,W2,Co}
  minBICsM1 <- apply(BICList[[i]][,-3],MARGIN=1,FUN=function(x) which.min(x))
  minAICsM1 <- apply(AICList[[i]][,-3],MARGIN=1,FUN=function(x) which.min(x))
  M1SelectBIC[i] <- round(sum(minBICsM1==3)/length(minBICsM1),2)
  M1SelectAIC[i] <- sum(minAICsM1==3)/length(minAICsM1)
  # choose model based on BIC & AIC alone for # M2={W1,W2,Co}
  minBICsM2 <- apply(BICList[[i]][,-4],MARGIN=1,FUN=function(x) which.min(x))
  minAICsM2 <- apply(AICList[[i]][,-4],MARGIN=1,FUN=function(x) which.min(x))
  M2SelectBIC[i] <- round(sum(minBICsM2==3)/length(minBICsM2),2)
  M2SelectAIC[i] <- sum(minAICsM2==3)/length(minAICsM2)
}
M1SelectBIC
M1SelectAIC
M2SelectBIC
M2SelectAIC

# now apply ECIC for all models for the first observed data set for each sample size
splinesW1ECIC <- list()
splinesW2ECIC <- list()
splinesClECIC <- list()
splinesCoECIC <- list()
modNamesW1 <- paste("splineW1,n=",ns,sep="")
modNamesW2 <- paste("splineW2,n=",ns,sep="")
modNamesCl <- paste("splineCl,n=",ns,sep="")
modNamesCo <- paste("splineCo,n=",ns,sep="")
# ecicModel() extracts info. from the input object and lists them 
# in a new ecicModel object
for(i in 1:length(ns))
{
  w1Temp <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfW1,degree=3))
  w2Temp <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfW2,degree=3))
  clTemp <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfCl,degree=3))
  coTemp <- lm(datList[[i]][,1]~bs(xs[[i]],df=dfCo,degree=3))
  splinesW1ECIC[[i]] <- ecicModel(model.name=w1Temp,ID=modNamesW1[i])
  splinesW2ECIC[[i]] <- ecicModel(model.name=w2Temp,ID=modNamesW2[i])
  splinesClECIC[[i]] <- ecicModel(model.name=clTemp,ID=modNamesCl[i])
  splinesCoECIC[[i]] <- ecicModel(model.name=coTemp,ID=modNamesCo[i])
}

# create model sets by sample size for M1 and M2
# ecicModelList() essentially throws errors if an invalid
# model set is created
modelSetsM1 <- list()
modelSetsM2 <- list()
for(i in 1:length(ns))
{
  modelSetsM1[[i]] <- ecicModelList(list(splinesW1ECIC[[i]],splinesW2ECIC[[i]],splinesCoECIC[[i]]))
  modelSetsM2[[i]] <- ecicModelList(list(splinesW1ECIC[[i]],splinesW2ECIC[[i]],splinesClECIC[[i]]))
}

#note that ecicModelList() is invoked in the ECIC function so the 
# above step isn't necessary 
# apply ECIC to spline models to see what happens asymptotically
# data=response variable?
# methods(IC) only shows IC.AIC, no BIC programmed?
# many warnings at the end
ECICResultsM1AIC <- list()
ECICResultsM2AIC <- list()
for(i in 1:length(ns))
{
  ECICResultsM1AIC[[i]] <- ECIC1.1Dev::ECIC(models=modelSetsM1[[i]],data=datList[[i]][,1],
                                 alpha=c(0.01,.05,.1),N=1000,ic="AIC")
  ECICResultsM2AIC[[i]] <- ECIC1.1Dev::ECIC(models=modelSetsM2[[i]],data=datList[[i]][,1],
                                      alpha=c(0.01,.05,.1),N=1000,ic="AIC")
  print(i)
}

# print decisions
for(i in 1:length(ns))
{
  print(ECICResultsM1AIC[[i]]$decisions)
  print(ECICResultsM2AIC[[i]]$decisions)
}