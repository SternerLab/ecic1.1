#####################################################################
#                       Description
#####################################################################

# This script assesses asymptotic behavior of ECIC under correct and 
# incorrect model specification
# for a model set of three Bernoulli distributions with parameters
# 0.5, 0.65, and 0.75
# We simulate 2 scenarios:
# Scenario #1: The true probability distribution is a Bernoulli 
# distribution with p=0.65 (correct model specification)
# Scenario #2: The true probability distribution is a Bernoulli
# distribution with p=0.68 (incorrect model specification)
# Note that we will use a simplified version of ECIC since we will
# take the parameters for each model to be fixed and
# thus do not need to undergo MLE estimation and correction as in the
# original algorithm
# Note that we will use the negative LL as the IC because of the fixed parameters

#####################################################################
#                       Functions
#####################################################################

# compute the log likelihood for a Bernoulli distribution 
LLBern <- function(dat,p)
{
  LL <- log(p)*sum(dat) + log(1-p)*sum(1-dat)
  return(LL)
}

# compute the BIC under a Bernoulli distribution
BICBern <- function(dat,p)
{
  LL <- log(p)*sum(dat) + log(1-p)*sum(1-dat)
  BIC <- -2*LL + 1*log(length(dat))
  return(BIC)
}

# compute the AIC under a Bernoulli distribution
AICBern <- function(dat,p)
{
  LL <- log(p)*sum(dat) + log(1-p)*sum(1-dat)
  AIC <- -2*LL + 2*1
  return(AIC)
}

# compute the difference in goodness of fit statistic for the DGOF simulation
# in step 4c of ECIC
DGOFSimComp <- function(ICScores,MbInd)
{
  DGOF <- ICScores[MbInd] - min(ICScores[-MbInd])
  return(unname(DGOF))
}


#####################################################################
#                      Load Packages/Data Generation
#####################################################################

# store parameters for the model set 
M <- c(0.5,0.65,0.75)
# Create 2 lists of vectors that hold data generated from
# a Bernoulli distribution
# the first list will possess data generated from a Bernoulli distribution
# with p=0.65 and the second with p=0.68
# Each element of each list will correspond to different sample sizes (see ns below)
# store the true parameter for scenario #1 and scenario #2
truePS1 <- 0.65
truePS2 <- 0.68
# set different sample sizes
ns <- c(20,50,100,200,500,1000)
datList1 <- list()
datList2 <- list()
set.seed(222)
for(i in 1:length(ns))
{
  tempN <- ns[i]
  # simulate data for Scenario #1
  # draw from a Bernoulli dist. with current n and p=truePS1
  tempDraw <- rbinom(n=tempN,size=1,prob=truePS1)
  # store vector of data of size n
  datList1[[i]] <- tempDraw
  # simulate data for Scenario #2
  # draw from a Bernoulli dist with current n & p=truePS2
  tempDraw <- rbinom(n=tempN,size=1,prob=truePS2)
  # store vector of data of size n
  datList2[[i]] <- tempDraw
}
# set names for datLists
names(datList1) <- paste("simulated w/ p=0.65,n=",ns,sep="")
names(datList2) <- paste("simulated w/ p=0.68,n=",ns,sep="")

#####################################################################
#                      Run Study for Model Set 1
#####################################################################

# ECIC step #1
# compute the IC under each model in the model set for each sample size
ICComps <- list()
for(i in 1:length(ns)) # i indexes sample size, j indexes elements in the model set
{
  tempICs <- rep(NA,length(M)) # vector w/ same no. of elements as model set to store -LL for each model
  for(j in 1:length(M))
  {
    tempICs[j] <- -1*LLBern(datList1[[i]],p=M[j]) 
  }
  names(tempICs) <- paste("LL,p=",M,sep="")
  ICComps[[i]] <- tempICs
}
names(ICComps) <- paste("true p=0.65,n=",ns,sep="")

# ECIC step #2
# apply which.min(x) to each element of the ICComps list to determine which
# element of the model set is observed as best (minimum negative LL) for each sample size
# note 1 corresponds to p=0.5, 2 to p=0.65, and 3 to 0.75
Mb <- unlist(lapply(ICComps,FUN=function(x) which.min(x)))
Mb

# use to store the observed DGOF for each sample size
obsDGOFs <- rep(NA,length(ns))
# use to store the alternative models i.e. model set without the best observed model
altMList <- list()
# ECIC step #3
# compute the observed DGOF and store the alternative models
for(i in 1:length(ns))
{
  # find index with lowest IC score
  tempBestInd <- Mb[i]
  # store remaining IC scores when the lowest one is removed
  tempICBestRemoved <- ICComps[[i]][-tempBestInd]
  # compute and store the observed DGOF
  obsDGOFs[i] <- unname(ICComps[[i]][tempBestInd]-min(tempICBestRemoved))
  # store the alternative model set
  altMList[[i]] <- M[-tempBestInd]
}

# sample size for estimating the probability of choosing the observed best model under the assumption an
# alternative model is true
N1 <- 2000
# sample size for simulating the DGOF distribution under the assumption that an alternative model is true
N2 <- 4000
# pre-specified type-1 error rate
alpha <- 0.05
# vector to store the minimum quantiles for each sample size from the simulated DGOF distributions
minQuantiles <- rep(NA,length(ns))
# list to store the simulated DGOF distributions for each sample size
DGOFList <- list()
# matrix to store the estimated probability that Mb will be observed best if Mi is assumed to be true
piHatiMat <- matrix(NA,nrow=length(ns),ncol=length(M[-1]))
set.seed(223)
# ECIC steps #4 and #5
# go through each of the alternative models, assume they are true, and compute quantiles
# i indexes sample size, j indexes models in the alternative set, and k indexes models in the full set 
for(i in 1:length(ns))
{
  # current sample size
  tempN <- ns[i]
  # alternate models for current sample size
  altModelsTemp <- altMList[[i]]
  # vector to store the estimated quantiles from the simulated DGOF distributions
  quantileVec <- rep(NA,length(altModelsTemp))
  # store the -LL for the best observed model for sample size ns[i]
  fMbTemp <- ICComps[[i]][Mb[i]]
  # store the index for the best observed model for sample size ns[i]
  fMbIndTemp <- Mb[i]
  # matrix to store the simulated DGOF distribution for each sample size
  # rows correspond to draws 1,...,N2 and cols correspond to the 2 alternative models
  DGOFDistsTemp <- matrix(NA,nrow=N2,ncol=length(altModelsTemp)) 
  for(j in 1:length(altModelsTemp))
  {
    # set the alternative model assumed to be true for this iteration
    tempTrueP <- altModelsTemp[j]
    # sample N1 new datasets of size tempN under the assumption that tempTrueP is the true parameter
    tempSimDraws1 <- rbinom(n=tempN*N1,size=1,prob=tempTrueP)
    # store the draws as a matrix
    tempSimDraws1 <- matrix(data=tempSimDraws1,nrow=tempN,ncol=N1)
    # compute the frequency at which the observed best model is selected as best under the assumption that
    # an alternative model is true
    # tempMat1 will be a # 3XN1 matrix where the rows corresponds to the models in the full model set
    # and the columns correspond to the simulated draws assuming tempTrueP is the true parameter
    tempMat1 <- NULL
    # iterate through each model in the full set and compute the IC for all N1 simulated draws
    for(k in 1:length(M)) 
    {
      tempICs <- -1*apply(X=tempSimDraws1,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempMat1 <- rbind(tempMat1,tempICs)
    }
    rownames(tempMat1) <- paste("ICVal,p=",M,sep="")
    # determine the model with the minimum IC for each draw
    tempBestMods <- apply(tempMat1,MARGIN=2,FUN=function(x) which.min(x))
    # determine relative frequency that the observed best model is selected as best in this simulated scenario
    tempPiHat <- sum(tempBestMods==fMbIndTemp)/length(tempBestMods)
    # store this probability in a matrix to look at them later
    piHatiMat[i,j] <- tempPiHat
    # now simulate the DGOF distribution by drawing N2 new datasets of size tempN under the assumption that 
    # tempTrueP is the true parameter
    tempSimDraws2 <- rbinom(n=tempN*N2,size=1,prob=tempTrueP)
    tempSimDraws2 <- matrix(data=tempSimDraws2,nrow=tempN,ncol=N2)
    # tempMat2 will be a # 3XN2 matrix where the rows corresponds to the models in the full model set
    # and the columns correspond to the simulated draws assuming tempTrueP is the true parameter
    tempMat2 <- NULL
    # iterate through each model in the full set and compute the IC for all N2 simulated draws
    for(k in 1:length(M)) 
    {
      tempICs <- -1*apply(X=tempSimDraws2,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempMat2 <- rbind(tempMat2,tempICs)
    }
    rownames(tempMat2) <- paste("BICVal,p=",M,sep="")
    # estimate tau 
    tauHati <- min(c(alpha/tempPiHat,1))
    # simulate DGOF values 
    tempDGOFs <- apply(tempMat2,MARGIN=2,FUN=function(x) DGOFSimComp(x,MbInd=fMbIndTemp))
    # store simulated DGOF values to look at them later
    DGOFDistsTemp[,j] <- tempDGOFs
    # estimate quantile value for the simulated DGOF distribution
    tempQuantile <- quantile(tempDGOFs,probs=tauHati)
    # plot the DGOF distribution 
    hist(tempDGOFs,main=paste("Estimated DGOF Dist. (Red Line=Est. Quantile)\n","Assumed True p=",tempTrueP,", n=",tempN)
          )
    abline(v=tempQuantile,col="red")
    quantileVec[j] <- tempQuantile
  }
  # store the simulated DGOF distributions for both alternative models for the given sample size i
  DGOFList[[i]] <- DGOFDistsTemp
  # store the minimum quantile produced from simulating the DGOF distributions under the assumption of the 
  # alternative models
  minQuantiles[i] <- min(quantileVec)
  print(i)
}

# print the estimated minimum quantiles
minQuantiles
# print the observed DGOFs
obsDGOFs
# combine vectors to set the ylim parameter in ensuing plot
yLimVec <- c(minQuantiles,obsDGOFs)
plot(obsDGOFs,type="l",xaxt="n",xlab="n", ylim=c(min(yLimVec),max(yLimVec)))
title(main="Observed DGOFs (Black) and Decision Threshold (Red) (Top Axis=Obs. Best Model)",line=3)
axis(side=1,at=1:length(ns),labels=ns)
lines(minQuantiles,col="red")
axis(side=3,at=1:length(ns),labels=Mb)


#####################################################################
#                      Run Study for Model Set 2
#####################################################################


# ECIC step #1
# compute the IC under each model in the model set for each sample size
ICComps <- list()
for(i in 1:length(ns)) # i indexes sample size, j indexes elements in the model set
{
  tempICs <- rep(NA,length(M)) # vector w/ same no. of elements as model set to store -LL for each model
  for(j in 1:length(M))
  {
    tempICs[j] <- -1*LLBern(datList2[[i]],p=M[j]) 
  }
  names(tempICs) <- paste("LL,p=",M,sep="")
  ICComps[[i]] <- tempICs
}
names(ICComps) <- paste("true p=0.65,n=",ns,sep="")

# ECIC step #2
# apply which.min(x) to each element of the ICComps list to determine which
# element of the model set is observed as best (minimum negative LL) for each sample size
# note 1 corresponds to p=0.5, 2 to p=0.65, and 3 to 0.75
Mb <- unlist(lapply(ICComps,FUN=function(x) which.min(x)))
Mb

# use to store the observed DGOF for each sample size
obsDGOFs <- rep(NA,length(ns))
# use to store the alternative models i.e. model set without the best observed model
altMList <- list()
# ECIC step #3
# compute the observed DGOF and store the alternative models
for(i in 1:length(ns))
{
  # find index with lowest IC score
  tempBestInd <- Mb[i]
  # store remaining IC scores when the lowest one is removed
  tempICBestRemoved <- ICComps[[i]][-tempBestInd]
  # compute and store the observed DGOF
  obsDGOFs[i] <- unname(ICComps[[i]][tempBestInd]-min(tempICBestRemoved))
  # store the alternative model set
  altMList[[i]] <- M[-tempBestInd]
}

# sample size for estimating the probability of choosing the observed best model under the assumption an
# alternative model is true
N1 <- 2000
# sample size for simulating the DGOF distribution under the assumption that an alternative model is true
N2 <- 4000
# pre-specified type-1 error rate
alpha <- 0.05
# vector to store the minimum quantiles for each sample size from the simulated DGOF distributions
minQuantiles <- rep(NA,length(ns))
# list to store the simulated DGOF distributions for each sample size
DGOFList <- list()
# matrix to store the estimated probability that Mb will be observed best if Mi is assumed to be true
piHatiMat <- matrix(NA,nrow=length(ns),ncol=length(M[-1]))
set.seed(223)
# ECIC steps #4 and #5
# go through each of the alternative models, assume they are true, and compute quantiles
# i indexes sample size, j indexes models in the alternative set, and k indexes models in the full set 
for(i in 1:length(ns))
{
  # current sample size
  tempN <- ns[i]
  # alternate models for current sample size
  altModelsTemp <- altMList[[i]]
  # vector to store the estimated quantiles from the simulated DGOF distributions
  quantileVec <- rep(NA,length(altModelsTemp))
  # store the -LL for the best observed model for sample size ns[i]
  fMbTemp <- ICComps[[i]][Mb[i]]
  # store the index for the best observed model for sample size ns[i]
  fMbIndTemp <- Mb[i]
  # matrix to store the simulated DGOF distribution for each sample size
  # rows correspond to draws 1,...,N2 and cols correspond to the 2 alternative models
  DGOFDistsTemp <- matrix(NA,nrow=N2,ncol=length(altModelsTemp)) 
  for(j in 1:length(altModelsTemp))
  {
    # set the alternative model assumed to be true for this iteration
    tempTrueP <- altModelsTemp[j]
    # sample N1 new datasets of size tempN under the assumption that tempTrueP is the true parameter
    tempSimDraws1 <- rbinom(n=tempN*N1,size=1,prob=tempTrueP)
    # store the draws as a matrix
    tempSimDraws1 <- matrix(data=tempSimDraws1,nrow=tempN,ncol=N1)
    # compute the frequency at which the observed best model is selected as best under the assumption that
    # an alternative model is true
    # tempMat1 will be a # 3XN1 matrix where the rows corresponds to the models in the full model set
    # and the columns correspond to the simulated draws assuming tempTrueP is the true parameter
    tempMat1 <- NULL
    # iterate through each model in the full set and compute the IC for all N1 simulated draws
    for(k in 1:length(M)) 
    {
      tempICs <- -1*apply(X=tempSimDraws1,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempMat1 <- rbind(tempMat1,tempICs)
    }
    rownames(tempMat1) <- paste("ICVal,p=",M,sep="")
    # determine the model with the minimum IC for each draw
    tempBestMods <- apply(tempMat1,MARGIN=2,FUN=function(x) which.min(x))
    # determine relative frequency that the observed best model is selected as best in this simulated scenario
    tempPiHat <- sum(tempBestMods==fMbIndTemp)/length(tempBestMods)
    # store this probability in a matrix to look at them later
    piHatiMat[i,j] <- tempPiHat
    # now simulate the DGOF distribution by drawing N2 new datasets of size tempN under the assumption that 
    # tempTrueP is the true parameter
    tempSimDraws2 <- rbinom(n=tempN*N2,size=1,prob=tempTrueP)
    tempSimDraws2 <- matrix(data=tempSimDraws2,nrow=tempN,ncol=N2)
    # tempMat2 will be a # 3XN2 matrix where the rows corresponds to the models in the full model set
    # and the columns correspond to the simulated draws assuming tempTrueP is the true parameter
    tempMat2 <- NULL
    # iterate through each model in the full set and compute the IC for all N2 simulated draws
    for(k in 1:length(M)) 
    {
      tempICs <- -1*apply(X=tempSimDraws2,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempMat2 <- rbind(tempMat2,tempICs)
    }
    rownames(tempMat2) <- paste("BICVal,p=",M,sep="")
    # estimate tau 
    tauHati <- min(c(alpha/tempPiHat,1))
    # simulate DGOF values 
    tempDGOFs <- apply(tempMat2,MARGIN=2,FUN=function(x) DGOFSimComp(x,MbInd=fMbIndTemp))
    # store simulated DGOF values to look at them later
    DGOFDistsTemp[,j] <- tempDGOFs
    # estimate quantile value for the simulated DGOF distribution
    tempQuantile <- quantile(tempDGOFs,probs=tauHati)
    # plot the DGOF distribution 
    hist(tempDGOFs,main=paste("Estimated DGOF Dist. (Red Line=Est. Quantile)\n","Assumed True p=",tempTrueP,", n=",tempN)
    )
    abline(v=tempQuantile,col="red")
    quantileVec[j] <- tempQuantile
  }
  # store the simulated DGOF distributions for both alternative models for the given sample size i
  DGOFList[[i]] <- DGOFDistsTemp
  # store the minimum quantile produced from simulating the DGOF distributions under the assumption of the 
  # alternative models
  minQuantiles[i] <- min(quantileVec)
  print(i)
}

# print the estimated minimum quantiles
minQuantiles
# print the observed DGOFs
obsDGOFs
# combine vectors to set the ylim parameter in ensuing plot
yLimVec <- c(minQuantiles,obsDGOFs)
plot(obsDGOFs,type="l",xaxt="n",xlab="n", ylim=c(min(yLimVec),max(yLimVec)))
title(main="Observed DGOFs (Black) and Decision Threshold (Red) (Top Axis=Obs. Best Model)",line=3)
axis(side=1,at=1:length(ns),labels=ns)
lines(minQuantiles,col="red")
axis(side=3,at=1:length(ns),labels=Mb)
