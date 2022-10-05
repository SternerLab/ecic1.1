#####################################################################
#                       Description
#####################################################################

# This script assesses asymptotic error probabilities for ECIC
# for a model set of three Bernoulli distributions with parameters
# 0.48,0.65, and 0.75
# We simulate 3 scenarios:
# Scenario #1: The true probability distribution is a Bernoulli 
# distribution with p=0.65
# Scenario #2: The true probability distribution is a Bernoulli
# distribution with p=0.68
# The purpose of this simulation study is to hopefully demonstrate that
# both the type 1 error rate (false positive rate) remains less than 
# or equal to a pre-specified value for all selected sample sizes and
# asymptotically approaches zero for both scenario #1 (true model 
# exists in the model set) and scenario #2 (true model does not 
# exist in the model set)
# Note that we will use a simplified version of ECIC since we will
# take the parameters for each model to be fixed and
# thus do not need to undergo MLE estimation and correction as in the
# original algorithm

#####################################################################
#                       Functions
#####################################################################

# compute the log likelihood for a Bernoulli distribution 
LLBern <- function(dat,p)
{
  LL <- log(p)*sum(dat) + log(1-p)*sum(1-dat)
  return(LL)
}

# a general function for computing the DGOF
DGOFGenComp <- function(ICScores)
{
  DGOF <- min(ICScores)-min(ICScores[-which.min(ICScores)])
  return(DGOF)
}

# compute the difference in goodness of fit statistic for the DGOF simulation
# in step 4c of ECIC
DGOFSimComp <- function(ICScores,MbInd)
{
  DGOF <- ICScores[MbInd] - min(ICScores[-MbInd])
  return(DGOF)
}

#####################################################################
#                      Load Packages/Data Generation
#####################################################################

# store parameters for the model set 
M <- c(0.48,0.65,0.75)
# store the cardinality of the model set
MLen <- length(M)
# Create 2 lists of matrices that hold data generated from
# a Bernoulli distribution
# the first list will possess data generated from a Bernoulli distribution
# with p=0.65 and the second with p=0.68
# Each element of each list will correspond to different sample sizes (see ns below)
# store the true parameter for scenario #1 and scenario #2
truePS1 <- 0.65
truePS2 <- 0.68
# set different sample sizes
ns <- c(5,15,20,25,30,40,50,70,150)
# store the cardinality of the sample sizes
nsLen <- length(ns)
# set the number of draws for each sample size
noDraws <- 2000
datList1 <- list()
datList2 <- list()
set.seed(222)
for(i in 1:nsLen)
{
  tempN <- ns[i]
  # simulate for model set 1
  # draw from a Bernoulli dist with current n & p
  tempDraws <- rbinom(n=tempN*noDraws,size=1,prob=truePS1)
  # convert vector of draws into a matrix of multiple draws
  tempDraws <- matrix(data=tempDraws,nrow=tempN,ncol=noDraws)
  # store matrix as the draws that are each of sample size tempN
  datList1[[i]] <- tempDraws
  colnames(datList1[[i]]) <- paste("Draw",1:noDraws,sep="")
  # simulate for model set 2
  # draw from a Bernoulli dist with current n & p
  tempDraws <- rbinom(n=tempN*noDraws,size=1,prob=truePS2)
  # convert vector of draws into a matrix of multiple draws
  tempDraws <- matrix(data=tempDraws,nrow=tempN,ncol=noDraws)
  # store matrix as the draws that are each of sample size tempN 
  datList2[[i]] <- tempDraws
  colnames(datList2[[i]]) <- paste("Draw",1:noDraws,sep="")
}
# set names for datLists
names(datList1) <- paste("True p=0.65,n=",ns,sep="")
names(datList2) <- paste("True p=0.65,n=",ns,sep="")

#####################################################################
#                      Run Study for Scenario 1
#####################################################################

# ECIC step #1
# compute the IC under each model in the model set for each sample size
ICComps <- list()
for(i in 1:nsLen)
{
  # MLenXnoDraws matrix holding the IC value for the 3 models for each draw
  tempMat <- matrix(NA,nrow=MLen,ncol=noDraws)
  for(j in 1:MLen)
  {
    tempICs <- -1*apply(X=datList1[[i]],MARGIN=2,FUN=function(x) LLBern(x,p=M[j]))
    tempMat[j,] <- tempICs
  }
  ICComps[[i]] <- tempMat
  # columns are draws and rows are the IC evaluated under different parameter values
  rownames(ICComps[[i]]) <- paste("IC Under p=",M,sep="")
  colnames(ICComps[[i]]) <- paste("Draw",1:noDraws,sep="")
}
names(ICComps) <- paste("True p=0.65,Draws of n=",ns,sep="")

# ECIC step #2
# apply which.min(x) across each column of each matrix in the ICComps list 
# to determine the observed best models
# list of length nsLen of 3XnoDraws vectors in each element that holds the lowest 
# IC score for each draw of each sample size
MbList <- list()
for(i in 1:nsLen)
{
  MbList[[i]] <- apply(X=ICComps[[i]],MARGIN=2,FUN=function(x) M[which.min(x)])
}
names(MbList) <- paste("Draws of n=",ns,",Mbs",sep="")

# make barplots for the percentage of the time each model was chosen as best for each sample size
for(i in 1:nsLen)
{
  # get the frequencies of models selected as best
  tempTable <- table(MbList[[i]])
  barplot(tempTable/noDraws,main=paste("Dist. of of Best Observed Models  (p=0.65 is true model)\n n=",ns[i]))
}

# list to store the observed best DGOF for each draw of each sample size
obsDGOFs <- list() 
# ECIC step #3
# compute the observed DGOFs
for(i in 1:nsLen)
{
  # compute the observed DGOFs for each draw of each sample size
  obsDGOFs[[i]] <- apply(X=ICComps[[i]],MARGIN=2,FUN=function(x) DGOFGenComp(x))
}
names(obsDGOFs) <- paste("Draws of n=",ns,",Observed DGOFs",sep="")

# sample size for estimating the probability of choosing the observed best model under the assumption an
# alternative model is true
N1 <- 2000
# sample size for simulating the DGOF distribution under the assumption that an alternative model is true
N2 <- 4000
# pre-specified type 1 error rate
alpha <- 0.05

# compute data from simulated models in ECIC step #4 in 
# advance to ease computational burden
# simDat1List holds the datasets to estimate \hat{\pi_i}=P_i(g(F)=M_b)
# simDat2List holds the datasets to estimate the DGOF distributions \Delta f_i
simDat1List <- list()
simDat2List <- list()
# initialize the above lists' elements as lists
for(i in 1:nsLen)
{
  simDat1List[[i]] <- list()
  simDat2List[[i]] <- list()
}

set.seed(19)
# simulate data for each sample size and probability in the model set
for(i in 1:nsLen) # i indexes sample sizes,
{
  tempN <- ns[i]
  for(j in 1:MLen) #  j indexes the assumed true probability
  {
    simDat1List[[i]][[j]] <- rbinom(n=tempN*N1,size=1,prob=M[j])
    simDat1List[[i]][[j]] <- matrix(data=simDat1List[[i]][[j]],nrow=tempN,ncol=N1)
    colnames(simDat1List[[i]][[j]]) <- paste("Draw",1:N1,sep="")
    simDat2List[[i]][[j]] <- rbinom(n=tempN*N2,size=1,prob=M[j])
    simDat2List[[i]][[j]] <- matrix(data=simDat2List[[i]][[j]],nrow=tempN,ncol=N2)
    colnames(simDat2List[[i]][[j]]) <- paste("Draw",1:N2,sep="")
  }
  names(simDat1List[[i]]) <- paste("True p=",M,sep="")
  names(simDat2List[[i]]) <- paste("True p=",M,sep="")
}
names(simDat1List) <- paste("Draws of n=",ns)
names(simDat2List) <- paste("Draws of n=",ns)

# each element in these lists will hold matrices of the -LL's under the 
# probabilities in the model sets for all draws of each sample size and 
# generating probability
ICsSimDat1 <- list()
ICsSimDat2 <- list()
# initialize the above lists' elements as lists 
for(i in 1:nsLen)
{
  ICsSimDat1[[i]] <- list()
  ICsSimDat2[[i]] <- list()
}

for(i in 1:nsLen) # i indexes sample sizes
{
  tempN <- ns[i]
  for(j in 1:MLen) # j indexes the assumed true probability 
  {
    # MLenXN1 matrix to store -LL's 
    tempMat1 <- matrix(NA,nrow=MLen,ncol=N1)
    colnames(tempMat1) <- paste("Draw",1:N1,sep="")
    # MLenXN2 matrix to store -LL's 
    tempMat2 <- matrix(NA,nrow=MLen,ncol=N2)
    colnames(tempMat2) <- paste("Draw",1:N2,sep="")
    tempDraws1 <- simDat1List[[i]][[j]] 
    tempDraws2 <- simDat2List[[i]][[j]] 
    for(k in 1:MLen) # compute -LL's under different parameters in the model set
    {
      tempICs1 <- -1*apply(X=tempDraws1,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempICs2 <- -1*apply(X=tempDraws2,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempMat1[k,] <- tempICs1
      tempMat2[k,] <- tempICs2
    }
    ICsSimDat1[[i]][[j]] <- tempMat1
    ICsSimDat2[[i]][[j]] <- tempMat2
    rownames(ICsSimDat1[[i]][[j]]) <- paste("-LL Under p=",M,sep="")
    rownames(ICsSimDat2[[i]][[j]]) <- paste("-LL Under p=",M,sep="")
  }
  names(ICsSimDat1[[i]]) <- paste("True p=",M,sep="")
  names(ICsSimDat2[[i]]) <- paste("True p=",M,sep="")
}
names(ICsSimDat1) <- paste("Draws of n=",ns,sep="")
names(ICsSimDat2) <- paste("Draws of n=",ns,sep="")

# determine the models with the minimum IC for each set of draws
minICList <- list()
for(i in 1:nsLen)
{
  minICList[[i]] <- list()
}

for(i in 1:nsLen) # i indexes sample size
{
  for(j in 1:MLen) # j indexes assumed true parameter
  {
    minICList[[i]][[j]] <- apply(ICsSimDat1[[i]][[j]],MARGIN=2,FUN=function(x) M[which.min(x)])
  }
  names(minICList[[i]]) <- paste("True p=",M,sep="")
}
names(minICList) <- paste("Draws of n=",ns,sep="")

# create a list of matrices that hold P_i(g(F)=M_b)
piHatList <- list()
for(i in 1:nsLen) # indexes the sample size
{
  tempMat <- matrix(data=NA,nrow=MLen,ncol=MLen)
  rownames(tempMat) <- paste("true p=",M,sep="")
  colnames(tempMat) <- paste("% of time p=",M," observed best ",sep="")
  for(j in 1:MLen) # indexes the true generating probability
  {
    for(k in 1:MLen) # indexes the times that the kth prob is selected as best
    {
      tempMat[j,k] <- sum(minICList[[i]][[j]]==M[k])/N1
    }
  }
  piHatList[[i]] <- tempMat
}
names(piHatList) <- paste("Draws of n=",ns,sep="")

# list to store the quantiles from alpha/P(Mb=M)
tauHatList <- list()
for(i in 1:nsLen)
{
  tauHatList[[i]] <- list()
}

for(i in 1:nsLen)
{
  tempMat <- piHatList[[i]]
  tempMat <- alpha/tempMat
  # if alpha/pi_hat >1, adjust the value down to 1
  tempMat[tempMat>1] <- 1
  colnames(tempMat) <- paste("tauth percentile using (p=",M," observed best) ",sep="")
  tauHatList[[i]] <- tempMat
}
names(tauHatList) <- paste("Draws of n=",ns,sep="")

# compute the DGOF distributions as in step 4c) of ECIC
DGOFList <- list()
for(i in 1:nsLen)
{
  DGOFList[[i]] <- list()
}

for(i in 1:nsLen) # index the sample size
{
  for(j in 1:MLen) # index the true generating probability 
  {
    tempDat <- ICsSimDat2[[i]][[j]]
    # matrix to store the DGOF distributions under the assumption that model k is
    # the observed best index
    tempMat <- matrix(NA,nrow=MLen,ncol=N2)
    colnames(tempMat) <- paste("Draw",1:N2,sep="")
    for(k in 1:MLen) #index the model assumed to be best
    {
      tempDGOFs <- apply(tempDat,MARGIN=2,FUN=function(x) DGOFSimComp(x,MbInd=k))
      tempMat[k,] <- tempDGOFs # store row as DGOF distribution
    }
    DGOFList[[i]][[j]] <- tempMat
    rownames(DGOFList[[i]][[j]]) <- paste("DGOFs,M_b=",M,sep="")
  }
  names(DGOFList[[i]]) <- paste("True p=",M,sep="")
}
names(DGOFList) <- paste("Draws of n=",ns,sep="")

# ECIC steps #4 and #5
# go through each of the alternative models, assume they are true, and compute quantiles
# i indexes sample size, j indexes models in the alternative set, and k indexes models in the full set
# list to store accept or reject decisions by sample size
aorRList <- list()
# list to store the decision thresholds
thresholds <- list()

for(i in 1:nsLen) # index sample size
{
  tempAorRVec <- rep(NA,noDraws)
  tempMatrix <- matrix(NA,nrow=noDraws,ncol=MLen-1) 
  for(j in 1:noDraws) #go through each draw for each sample size and perform ECIC
  {
    # identify the observed best model
    tempMbInd <- which(M==MbList[[i]][j])
    # identify the observed DGOF
    tempObsDGOF <- obsDGOFs[[i]][j]
    # identify the alternative models
    altModels <- M[-tempMbInd]
    tempDGOFQuantiles <- rep(NA,MLen-1)
    for(k in 1:(MLen-1)) # assume k is index for the true model
    {
      # now just retrieve the quantile and DGOF distribution
      # to compare it to the observed DGOF
      curAltModel <- altModels[k]
      curAltModelInd <- which(M==curAltModel)
      tempTau <- tauHatList[[i]][curAltModelInd,tempMbInd]
      tempDGOFs <- DGOFList[[i]][[curAltModelInd]][tempMbInd,]
      tempDGOFQuantiles[k] <- quantile(tempDGOFs,probs=tempTau)
    }
    tempMatrix[j,] <- tempDGOFQuantiles
    # take the alternative model quantile estimates
    tempFinQuantile <- min(tempDGOFQuantiles)
    # store 1 if observed model is chosen as best and 0 otherwise
    tempAorRVec[j] <- ifelse(test=tempObsDGOF<tempFinQuantile,yes=1,no=0)
  }
  thresholds[[i]] <- tempMatrix
  aorRList[[i]] <- tempAorRVec
  print(i)
}

# assess observed best model with ECIC choice
assessList <- list()
for(i in 1:nsLen)
{
  assessList[[i]] <- rbind(MbList[[i]],aorRList[[i]])
}

# plot the percentage of the time no decision is made by sample size
fracNoDec <- sapply(aorRList,FUN=function(x) sum(x==0)/length(x))
plot(x=ns,y=fracNoDec,main="Percentage of Time No Decision is Made by Sample Size",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

# Look at percentage of times the true model was observed as best but no decision was made
# subset assessList to only contain observations where with no decision i.e. 0
noDecAssessList <- lapply(assessList,FUN=function(x) x[,x[2,]==0])
rightObsUndec <- sapply(noDecAssessList,FUN=function(x) sum(x[1,]==truePS1)/noDraws)
# note that once a decision is always made (which occurs asymptotically)
# there will not be a point corresponding to the sample size on this following plot
plot(x=ns,y=rightObsUndec,main="Percentage of Time True Model Observed When No Decision Made",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)
# take a look at type 1 error rate
# subset assess list by only when a decision was made i.e. second row = 1
decAssessList <- lapply(assessList,FUN=function(x) x[,x[2,]==1])
# compute type 1 error rate
t1ErrorRates <- sapply(decAssessList,FUN=function(x) sum(x[1,]!=truePS1)/noDraws)
plot(x=ns,y=t1ErrorRates,main="Type 1 Error Rates by Sample Size",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

# Due to the discreteness of the DGOF distribution, we cannot precisely control the 
# error rates. Let's see what happens when we use the use the next smallest value
# as the quantile for a more conservative estimate
aorRList <- list()
# list to store the decision thresholds
thresholds <- list()

for(i in 1:nsLen) # index sample size
{
  tempAorRVec <- rep(NA,noDraws)
  tempMatrix <- matrix(NA,nrow=noDraws,ncol=MLen-1) 
  for(j in 1:noDraws) #go through each draw for each sample size and perform ECIC
  {
    # identify the observed best model
    tempMbInd <- which(M==MbList[[i]][j])
    # identify the observed DGOF
    tempObsDGOF <- obsDGOFs[[i]][j]
    # identify the alternative models
    altModels <- M[-tempMbInd]
    tempDGOFQuantiles <- rep(NA,MLen-1)
    for(k in 1:(MLen-1)) # assume k is index for the true model
    {
      # now just retrieve the quantile and DGOF distribution
      # to compare it to the observed DGOF
      curAltModel <- altModels[k]
      curAltModelInd <- which(M==curAltModel)
      tempTau <- tauHatList[[i]][curAltModelInd,tempMbInd]
      tempDGOFs <- DGOFList[[i]][[curAltModelInd]][tempMbInd,]
      uniqueDGOFs <- sort(unique(tempDGOFs))
      tempQuant <- quantile(tempDGOFs,probs=tempTau)
      tempQuantInd <- which(uniqueDGOFs==tempQuant)
      # adjust quantile to next smallest value if there exists one
      if(tempQuantInd>1)
        tempQuant <- uniqueDGOFs[tempQuantInd-1]
      tempDGOFQuantiles[k] <- tempQuant
    }
    tempMatrix[j,] <- tempDGOFQuantiles
    # take the alternative model quantile estimates
    tempFinQuantile <- min(tempDGOFQuantiles)
    # store 1 if observed model is chosen as best and 0 otherwise
    tempAorRVec[j] <- ifelse(test=tempObsDGOF<tempFinQuantile,yes=1,no=0)
  }
  thresholds[[i]] <- tempMatrix
  aorRList[[i]] <- tempAorRVec
  print(i)
}

# assess observed best model with ECIC choice
assessList <- list()
for(i in 1:nsLen)
{
  assessList[[i]] <- rbind(MbList[[i]],aorRList[[i]])
}
# take a look at type 1 error rate
# subset assess list by only when a decision was made i.e. second row = 1
decAssessList <- lapply(assessList,FUN=function(x) x[,x[2,]==1])
# compute type 1 error rate
t1ErrorRates <- sapply(decAssessList,FUN=function(x) sum(x[1,]!=truePS1)/noDraws)
plot(x=ns,y=t1ErrorRates,main="Type 1 Error Rates by Sample Size",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

#####################################################################
#                      Run Study for Scenario 2
#####################################################################

# ECIC step #1
# compute the IC under each model in the model set for each sample size
ICComps <- list()
for(i in 1:nsLen)
{
  # MLenXnoDraws matrix holding the IC value for the 3 models for each draw
  tempMat <- matrix(NA,nrow=MLen,ncol=noDraws)
  for(j in 1:MLen)
  {
    tempICs <- -1*apply(X=datList2[[i]],MARGIN=2,FUN=function(x) LLBern(x,p=M[j]))
    tempMat[j,] <- tempICs
  }
  ICComps[[i]] <- tempMat
  # columns are draws and rows are the IC evaluated under different parameter values
  rownames(ICComps[[i]]) <- paste("IC Under p=",M,sep="")
  colnames(ICComps[[i]]) <- paste("Draw",1:noDraws,sep="")
}
names(ICComps) <- paste("True p=0.65,Draws of n=",ns,sep="")

# ECIC step #2
# apply which.min(x) across each column of each matrix in the ICComps list 
# to determine the observed best models
# list of length nsLen of 3XnoDraws vectors in each element that holds the lowest 
# IC score for each draw of each sample size
MbList <- list()
for(i in 1:nsLen)
{
  MbList[[i]] <- apply(X=ICComps[[i]],MARGIN=2,FUN=function(x) M[which.min(x)])
}
names(MbList) <- paste("Draws of n=",ns,",Mbs",sep="")

# make barplots for the percentage of the time each model was chosen as best for each sample size
for(i in 1:nsLen)
{
  # get the frequencies of models selected as best
  tempTable <- table(MbList[[i]])
  barplot(tempTable/noDraws,main=paste("Dist. of of Best Observed Models  (p=0.65 is true model)\n n=",ns[i]))
}

# list to store the observed best DGOF for each draw of each sample size
obsDGOFs <- list() 
# ECIC step #3
# compute the observed DGOFs
for(i in 1:nsLen)
{
  # compute the observed DGOFs for each draw of each sample size
  obsDGOFs[[i]] <- apply(X=ICComps[[i]],MARGIN=2,FUN=function(x) DGOFGenComp(x))
}
names(obsDGOFs) <- paste("Draws of n=",ns,",Observed DGOFs",sep="")

# sample size for estimating the probability of choosing the observed best model under the assumption an
# alternative model is true
N1 <- 2000
# sample size for simulating the DGOF distribution under the assumption that an alternative model is true
N2 <- 4000
# pre-specified type 1 error rate
alpha <- 0.05

# compute data from simulated models in ECIC step #4 in 
# advance to ease computational burden
# simDat1List holds the datasets to estimate \hat{\pi_i}=P_i(g(F)=M_b)
# simDat2List holds the datasets to estimate the DGOF distributions \Delta f_i
simDat1List <- list()
simDat2List <- list()
# initialize the above lists' elements as lists
for(i in 1:nsLen)
{
  simDat1List[[i]] <- list()
  simDat2List[[i]] <- list()
}

set.seed(19)
# simulate data for each sample size and probability in the model set
for(i in 1:nsLen) # i indexes sample sizes,
{
  tempN <- ns[i]
  for(j in 1:MLen) #  j indexes the assumed true probability
  {
    simDat1List[[i]][[j]] <- rbinom(n=tempN*N1,size=1,prob=M[j])
    simDat1List[[i]][[j]] <- matrix(data=simDat1List[[i]][[j]],nrow=tempN,ncol=N1)
    colnames(simDat1List[[i]][[j]]) <- paste("Draw",1:N1,sep="")
    simDat2List[[i]][[j]] <- rbinom(n=tempN*N2,size=1,prob=M[j])
    simDat2List[[i]][[j]] <- matrix(data=simDat2List[[i]][[j]],nrow=tempN,ncol=N2)
    colnames(simDat2List[[i]][[j]]) <- paste("Draw",1:N2,sep="")
  }
  names(simDat1List[[i]]) <- names(simDat2List[[i]]) <- paste("True p=",M,sep="")
}
names(simDat1List) <- paste("Draws of n=",ns)
names(simDat2List) <- paste("Draws of n=",ns)

# each element in these lists will hold matrices of the -LL's under the 
# probabilities in the model sets for all draws of each sample size and 
# generating probability
ICsSimDat1 <- list()
ICsSimDat2 <- list()
# initialize the above lists' elements as lists 
for(i in 1:nsLen)
{
  ICsSimDat1[[i]] <- list()
  ICsSimDat2[[i]] <- list()
}

for(i in 1:nsLen) # i indexes sample sizes
{
  tempN <- ns[i]
  for(j in 1:MLen) # j indexes the assumed true probability 
  {
    # MLenXN1 matrix to store -LL's 
    tempMat1 <- matrix(NA,nrow=MLen,ncol=N1)
    colnames(tempMat1) <- paste("Draw",1:N1,sep="")
    # MLenXN2 matrix to store -LL's 
    tempMat2 <- matrix(NA,nrow=MLen,ncol=N2)
    colnames(tempMat2) <- paste("Draw",1:N2,sep="")
    tempDraws1 <- simDat1List[[i]][[j]] 
    tempDraws2 <- simDat2List[[i]][[j]] 
    for(k in 1:MLen) # compute -LL's under different parameters in the model set
    {
      tempICs1 <- -1*apply(X=tempDraws1,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempICs2 <- -1*apply(X=tempDraws2,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
      tempMat1[k,] <- tempICs1
      tempMat2[k,] <- tempICs2
    }
    ICsSimDat1[[i]][[j]] <- tempMat1
    ICsSimDat2[[i]][[j]] <- tempMat2
    rownames(ICsSimDat1[[i]][[j]]) <- paste("-LL Under p=",M,sep="")
    rownames(ICsSimDat2[[i]][[j]]) <- paste("-LL Under p=",M,sep="")
  }
  names(ICsSimDat1[[i]]) <- paste("True p=",M,sep="")
  names(ICsSimDat2[[i]]) <- paste("True p=",M,sep="")
}
names(ICsSimDat1) <- paste("Draws of n=",ns,sep="")
names(ICsSimDat2) <- paste("Draws of n=",ns,sep="")

# determine the models with the minimum IC for each set of draws
minICList <- list()
for(i in 1:nsLen)
{
  minICList[[i]] <- list()
}

for(i in 1:nsLen) # i indexes sample size
{
  for(j in 1:MLen) # j indexes assumed true parameter
  {
    minICList[[i]][[j]] <- apply(ICsSimDat1[[i]][[j]],MARGIN=2,FUN=function(x) M[which.min(x)])
  }
  names(minICList[[i]]) <- paste("True p=",M,sep="")
}
names(minICList) <- paste("Draws of n=",ns,sep="")

# create a list of matrices that hold P_i(g(F)=M_b)
piHatList <- list()
for(i in 1:nsLen) # indexes the sample size
{
  tempMat <- matrix(data=NA,nrow=MLen,ncol=MLen)
  rownames(tempMat) <- paste("true p=",M,sep="")
  colnames(tempMat) <- paste("% of time p=",M," observed best ",sep="")
  for(j in 1:MLen) # indexes the true generating probability
  {
    for(k in 1:MLen) # indexes the times that the kth prob is selected as best
    {
      tempMat[j,k] <- sum(minICList[[i]][[j]]==M[k])/N1
    }
  }
  piHatList[[i]] <- tempMat
}
names(piHatList) <- paste("Draws of n=",ns,sep="")

# list to store the quantiles from alpha/P(Mb=M)
tauHatList <- list()
for(i in 1:nsLen)
{
  tauHatList[[i]] <- list()
}

for(i in 1:nsLen)
{
  tempMat <- piHatList[[i]]
  tempMat <- alpha/tempMat
  # if alpha/pi_hat >1, adjust the value down to 1
  tempMat[tempMat>1] <- 1
  colnames(tempMat) <- paste("tauth percentile using (p=",M," observed best) ",sep="")
  tauHatList[[i]] <- tempMat
}
names(tauHatList) <- paste("Draws of n=",ns,sep="")

# compute the DGOF distributions as in step 4c) of ECIC
DGOFList <- list()
for(i in 1:nsLen)
{
  DGOFList[[i]] <- list()
}

for(i in 1:nsLen) # index the sample size
{
  for(j in 1:MLen) # index the true generating probability 
  {
    tempDat <- ICsSimDat2[[i]][[j]]
    # matrix to store the DGOF distributions under the assumption that model k is
    # the observed best index
    tempMat <- matrix(NA,nrow=MLen,ncol=N2)
    colnames(tempMat) <- paste("Draw",1:N2,sep="")
    for(k in 1:MLen) #index the model assumed to be best
    {
      tempDGOFs <- apply(tempDat,MARGIN=2,FUN=function(x) DGOFSimComp(x,MbInd=k))
      tempMat[k,] <- tempDGOFs # store row as DGOF distribution
    }
    DGOFList[[i]][[j]] <- tempMat
    rownames(DGOFList[[i]][[j]]) <- paste("DGOFs,M_b=",M,sep="")
  }
  names(DGOFList[[i]]) <- paste("True p=",M,sep="")
}
names(DGOFList) <- paste("Draws of n=",ns,sep="")

# ECIC steps #4 and #5
# go through each of the alternative models, assume they are true, and compute quantiles
# i indexes sample size, j indexes models in the alternative set, and k indexes models in the full set
# list to store accept or reject decisions by sample size
aorRList <- list()
# list to store the decision thresholds
thresholds <- list()

for(i in 1:nsLen) # index sample size
{
  tempAorRVec <- rep(NA,noDraws)
  tempMatrix <- matrix(NA,nrow=noDraws,ncol=MLen-1) 
  for(j in 1:noDraws) #go through each draw for each sample size and perform ECIC
  {
    # identify the observed best model
    tempMbInd <- which(M==MbList[[i]][j])
    # identify the observed DGOF
    tempObsDGOF <- obsDGOFs[[i]][j]
    # identify the alternative models
    altModels <- M[-tempMbInd]
    tempDGOFQuantiles <- rep(NA,MLen-1)
    for(k in 1:(MLen-1)) # assume k is index for the true model
    {
      # now just retrieve the quantile and DGOF distribution
      # to compare it to the observed DGOF
      curAltModel <- altModels[k]
      curAltModelInd <- which(M==curAltModel)
      tempTau <- tauHatList[[i]][curAltModelInd,tempMbInd]
      tempDGOFs <- DGOFList[[i]][[curAltModelInd]][tempMbInd,]
      tempDGOFQuantiles[k] <- quantile(tempDGOFs,probs=tempTau)
    }
    tempMatrix[j,] <- tempDGOFQuantiles
    # take the alternative model quantile estimates
    tempFinQuantile <- min(tempDGOFQuantiles)
    # store 1 if observed model is chosen as best and 0 otherwise
    tempAorRVec[j] <- ifelse(test=tempObsDGOF<tempFinQuantile,yes=1,no=0)
  }
  thresholds[[i]] <- tempMatrix
  aorRList[[i]] <- tempAorRVec
  print(i)
}

# assess observed best model with ECIC choice
assessList <- list()
for(i in 1:nsLen)
{
  assessList[[i]] <- rbind(MbList[[i]],aorRList[[i]])
}

# plot the percentage of the time no decision is made by sample size
fracNoDec <- sapply(aorRList,FUN=function(x) sum(x==0)/noDraws)
plot(x=ns,y=fracNoDec,main="Percentage of Time No Decision is Made by Sample Size",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

# Look at percentage of times the true model was observed as best but no decision was made
# subset assessList to only contain observations where with no decision i.e. 0
noDecAssessList <- lapply(assessList,FUN=function(x) x[,x[2,]==0])
rightObsUndec <- sapply(noDecAssessList,FUN=function(x) sum(x[1,]==truePS1)/noDraws)
# note that once a decision is always made (which occurs asymptotically)
# there will not be a point corresponding to the sample size on this following plot
plot(x=ns,y=rightObsUndec,main="Percentage of Time True Model Observed When No Decision Made",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)
# take a look at type 1 error rate
# subset assess list by only when a decision was made i.e. second row = 1
decAssessList <- lapply(assessList,FUN=function(x) x[,x[2,]==1])
# compute type 1 error rate
t1ErrorRates <- sapply(decAssessList,FUN=function(x) sum(x[1,]!=truePS1)/noDraws)
plot(x=ns,y=t1ErrorRates,main="Type 1 Error Rates by Sample Size",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

