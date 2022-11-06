#####################################################################
#                       Functions
#####################################################################

# LL=-n/2*log(2\pi)-n/2*log(\sigma^2)-1/(2\sigma^2)\sum(x_j-\mu)^2
# AIC = -2LL + 2p
# BIC = -2LL + plog(n)
# compute the BIC for a normal model under the 3 model possibilities
BICNorm <- function(dat,model)
{
  n <- length(dat)
  if(model=="N(0,1)")
    BIC <- -2*(-n/2*log(2*pi)-n/2*log(1)-1/(2*1)*sum((dat-0)^2))
  else if(model=="N(mu,1)")
  {
    muMLE <- mean(dat)
    BIC <- -2*(-n/2*log(2*pi)-n/2*log(1)-1/(2*1)*sum((dat-muMLE)^2)) + 1*log(n)
  }
  else
  {
    muMLE <- mean(dat)
    sig2MLE <- (1/n)*sum((dat-muMLE)^2)
    BIC <- -2*(-n/2*log(2*pi)-n/2*log(sig2MLE)-1/(2*sig2MLE)*sum((dat-muMLE)^2)) + 2*log(n)
  }
  return(BIC)
}

# compute the LL for a normal model under the 3 model possibilities
LLNorm <- function(dat,model)
{
  n <- length(dat)
  if(model=="N(0,1)")
  {
    LL <- -n/2*log(2*pi)-n/2*log(1)-1/(2*1)*sum((dat-0)^2)
  }
  else if(model=="N(mu,1)")
  {
    muMLE <- mean(dat)
    LL <- -n/2*log(2*pi)-n/2*log(1)-1/(2*1)*sum((dat-muMLE)^2)
  }
  else
  {
    muMLE <- mean(dat)
    sig2MLE <- 1/n*sum((dat-muMLE)^2)
    LL <- -n/2*log(2*pi)-n/2*log(sig2MLE)-1/(2*sig2MLE)*sum((dat-muMLE)^2)
  }
  return(LL)
}

# a general function for computing the DGOF
DGOFGenComp <- function(ICScores)
{
  DGOF <- min(ICScores)-min(ICScores[-which.min(ICScores)])
  return(DGOF)
}



#####################################################################
#                      Load Packages/Data Generation
#####################################################################

library(Rcpp)
library(RcppArmadillo)
# load C++ script that is used for simulation steps in ECIC
sourceCpp("C:\\Users\\Djizi\\OneDrive\\Desktop\\Sterner_Project\\ECIC_Code\\Norm_Sim_2_RCppCode.cpp")

trueMu <- 0.5
trueSig <- 1
trueModel <- "N(0.5,1)"
M <- c("N(0,1)","N(mu,1)","N(mu,sigma)")
# model that should be selected
closestMod <- "N(mu,1)"
MLen <- length(M)
# set different sample sizes
ns <-  c(3,10,50,100,200,350) 
# store the cardinality of the sample sizes
nsLen <- length(ns)
# set the number of draws for each sample size
noDraws <- 300
datList <- list()
set.seed(225)
for(i in 1:nsLen)
{
  tempN <- ns[i]
  # simulate for model set 1
  # draw from a normal dist with current n & p
  tempDraws <- rnorm(n=tempN*noDraws,mean=trueMu,sd=trueSig)
  # convert vector of draws into a matrix of multiple draws
  tempDraws <- matrix(data=tempDraws,nrow=tempN,ncol=noDraws)
  datList[[i]] <- tempDraws
  colnames(datList[[i]]) <- paste("Draw",1:noDraws,sep="")
}
names(datList) <- paste("Draws of n=",ns,sep="")

# compute MLE estimates in advance
MLEList <- list()
for(i in 1:nsLen)
{
  MLEList[[i]] <- list()
}

for(i in 1:nsLen)
{
  tempN <- ns[i]
  tempDatList <- datList[[i]]
  tempMeanMLEs <- apply(X=tempDatList,MARGIN=2,FUN=function(x) mean(x))
  # use computational form of 1/n\sum(x_j-\overline{x})^2 i.e. \sum x_j^2/n - \overline{x}^2
  tempSig2MLEs <- apply(X=tempDatList,MARGIN=2,FUN=function(x) sum(x^2))/tempN - tempMeanMLEs^2
  # MLE's for N(0,1) (just 0,1 b/c no parameter estimation)
  MLEList[[i]][[1]] <- matrix(rep(c(0,1),noDraws),nrow=2,ncol=noDraws)
  rownames(MLEList[[i]][[1]]) <- c("fixedMean","fixedsig2")
  # MLE's for N(mu,1) 
  MLEList[[i]][[2]] <- rbind(tempMeanMLEs,1)
  rownames(MLEList[[i]][[2]]) <- c("MLEmean","fixedsig2")
  # MLE's for N(mu,sigma)
  MLEList[[i]][[3]] <- rbind(tempMeanMLEs,tempSig2MLEs)
  rownames(MLEList[[i]][[3]]) <- c("MLEmean","MLEsig2")
  colnames(MLEList[[i]][[1]]) <- colnames(MLEList[[i]][[2]]) <- colnames(MLEList[[i]][[3]]) <- paste("Draw",1:noDraws,sep="")
  names(MLEList[[i]]) <- paste("MLEs for model",M)
}
names(MLEList) <- paste("Draws of n=",ns,sep="")

#####################################################################
#             Begin ECIC
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
    tempICs <- apply(X=datList[[i]],MARGIN=2,FUN=function(x) BICNorm(x,M[j]))
    tempMat[j,] <- tempICs
  }
  ICComps[[i]] <- tempMat
  # columns are draws and rows are the IC evaluated under different parameter values
  rownames(ICComps[[i]]) <- paste("IC Under ",M,sep="")
  colnames(ICComps[[i]]) <- paste("Draw",1:noDraws,sep="")
}
names(ICComps) <- paste("True model:",trueModel,",Draws of n=",ns,sep="")

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
  barplot(tempTable/noDraws,main=paste("Dist. of Best Observed Models(mu=",trueMu,"sigma=",trueSig,"is true model)\n n=",ns[i]))
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
N1 <- 300
# sample size for simulating the DGOF distribution under the assumption that an alternative model is true
N2 <- 700
# pre-specified type-1 error rate
alpha <- 0.15

# compute data from simulated models in ECIC step #4 in 
# advance to ease computational burden
# simDat1List holds the datasets to estimate \hat{\pi_i}=P_i(g(F)=M_b)
# simDat2List holds the datasets to estimate the DGOF distributions \Delta f_i

set.seed(19)
ptm=proc.time() #Start timing
simDat1List <- normDatSim(ns,MLEList,N1)
simDat2List <- normDatSim(ns,MLEList,N2)
proc.time() - ptm
names(simDat1List) <- paste("Draws of n=",ns,sep="")
names(simDat2List) <- paste("Draws of n=",ns,sep="")
# provide names for both simulated sets
for(i in 1:nsLen)
{
  for(j in 1:MLen) # indexes the models in the model set
  {
    names(simDat1List[[i]][[j]]) <- paste("Obs",1:noDraws,";",N1,"Simulated Draws")
    names(simDat2List[[i]][[j]]) <- paste("Obs",1:noDraws,";",N2,"Simulated Draws")
  }
  names(simDat1List[[i]]) <- paste("Generated from",M)
  names(simDat2List[[i]]) <- paste("Generated from",M)
}



# Use C++ code to simulate datasets
# the RcpplmComps inputs the list of simulated draws from the MLE fits 
# and the model set. Matrices of BICs are returned 
ptm=proc.time() #Start timing 
ICsSimDat1 <- RcppICComps(simDat1List,M)
proc.time() - ptm
# label the elements in ICsSimDat1
for(i in 1:nsLen)
{
  for(j in 1:MLen)
  {
    for(k in 1:noDraws)
    {
      colnames(ICsSimDat1[[i]][[j]][[k]]) <- paste("BIC Under", M)
    }
    names(ICsSimDat1[[i]][[j]]) <- paste("lm Fit for Obs",1:noDraws)
  }
  names(ICsSimDat1[[i]]) <- paste("Generated from",M)
}
names(ICsSimDat1) <- paste("Draws of n=",ns,sep="")

# determine the model with the minimum IC for each set of draws
minICList <- list()
for(i in 1:nsLen)
{
  minICList[[i]] <- list()
  for(j in 1:MLen)
    minICList[[i]][[j]] <- list()
}

for(i in 1:nsLen) # i indexes sample size
{
  for(j in 1:MLen) # j indexes assumed true parameter
  {
    for(k in 1:noDraws)
    {
      minICList[[i]][[j]][[k]] <- apply(ICsSimDat1[[i]][[j]][[k]],MARGIN=1,FUN=function(x) M[which.min(x)])
    }
    names(minICList[[i]][[j]]) <- paste("Normal Fit for Obs",1:noDraws)
  }
  names(minICList[[i]]) <- paste("Generated from",M)
  print(i)
}
names(minICList) <- paste("Draws of n=",ns,sep="")

# create a list of matrices that hold P_i(g(F)=M_b)
piHatList <- list()
for(i in 1:nsLen)
{
  piHatList[[i]] <- list()
  for(k in 1:noDraws)
  {
    piHatList[[i]][[k]] <- list()
  }
  names(piHatList[[i]]) <- paste("Normal Fit for Obs",1:noDraws)
}
names(piHatList) <- names(piHatList) <- paste("n=",ns,sep="")

for(i in 1:nsLen) # indexes the sample size
{
  for(k in 1:noDraws) # indexes the times that the kth prob is selected as best
  {
    tempMat <- matrix(data=NA,nrow=MLen,ncol=MLen)
    rownames(tempMat) <- paste("true model=",M,sep="")
    colnames(tempMat) <- paste("% of time ",M," chosen ",sep="")
    for(d in 1:MLen) # index the true generating model
      for(f in 1:MLen) # index the model being considered 
        tempMat[d,f] <- sum(minICList[[i]][[d]][[k]]==M[f])/N1
    piHatList[[i]][[k]] <- tempMat
  }
  print(i)
}
names(piHatList) <- paste("n=",ns,sep="")

rm(ICsSimDat1)
rm(minICList)
gc()

# Use C++ code to simulate datasets
# the RcppSimDGOFs inputs the list of simulated draws from the MLE fits 
# and the model set. DGOF distributions are returned
DGOFList <- RcppSimDGOFs(simDat2List,M)
# label the elements in DGOFList
for(i in 1:nsLen)
{
  for(j in 1:MLen)
  {
    for(k in 1:noDraws)
    {
      colnames(DGOFList[[i]][[j]][[k]]) <- paste("DGOF Under", M," Observed Best")
    }
    names(DGOFList[[i]][[j]]) <- paste("Normal Fit for Obs",1:noDraws)
  }
  names(DGOFList[[i]]) <- paste("Generated from",M)
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
  for(k in 1:noDraws) # go through each draw for each sample size and perform ECIC
  {
    # identify the observed best model
    tempMbInd <- which(M==MbList[[i]][k])
    # identify the observed DGOF
    tempObsDGOF <- obsDGOFs[[i]][k]
    # identify the alternative models
    altModels <- M[-tempMbInd]
    tempDGOFQuantiles <- rep(NA,MLen-1)
    for(h in 1:(MLen-1)) # assume h is index for the true model
    {
      # now just retrieve the quantile and DGOF distribution
      # to compare it to the observed DGOF
      curAltModel <- altModels[h]
      curAltModelInd <- which(M==curAltModel)
      probThisFalseMbInd <- piHatList[[i]][[k]][curAltModelInd,tempMbInd]
      # this was prompted by issue where probThisFalseMbInd=0 and totalProbFalseMb=0
      # check with Beckett later if this is the right move...
      if(probThisFalseMbInd==0) tempTau <- 1
      else
      {
        totalProbFalseMb <- sum(piHatList[[i]][[k]][curAltModelInd,-curAltModelInd])
        tempTau <- alpha*probThisFalseMbInd/totalProbFalseMb
        if(tempTau>1) tempTau <- 1
      }
      tempDGOFs <- DGOFList[[i]][[curAltModelInd]][[k]][,tempMbInd]
      tempDGOFQuantiles[h] <- quantile(tempDGOFs,probs=tempTau)
    }
    tempMatrix[k,] <- tempDGOFQuantiles
    # take the alternative model quantile estimates
    tempFinQuantile <- min(tempDGOFQuantiles)
    # store 1 if observed model is chosen as best and 0 otherwise
    tempAorRVec[k] <- ifelse(test=tempObsDGOF<=tempFinQuantile,yes=1,no=0)
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

DecisionRates <- sapply(X=assessList,FUN=function(x) sum(as.numeric(x[2,]))/noDraws)
plot(x=ns,y=DecisionRates,main="Proportion of runs a model was selected",
     xlab="Sample Size",ylab="% of a model was selected",pch=16,ylim = c(0,1))
# take a look at type 1 error rate
# subset assess list by only when a decision was made i.e. second row = 1
decAssessList <- lapply(X=assessList,FUN=function(x) x[,x[2,]==1])
# compute type 1 error rate
t1ErrorRates <- sapply(X=decAssessList,FUN=function(x) sum(x[1,]!=closestMod)/noDraws)
plot(x=ns,y=t1ErrorRates,main=c("Type 1 Error Rates by Sample Size at alpha=",alpha),
     xlab="Sample Size",ylab="% of Time Wrong Best Model Selected",pch=16,ylim = c(0,2*alpha))
CorrectRates <- sapply(X=decAssessList,FUN=function(x) sum(x[1,]==closestMod)/noDraws)
plot(x=ns,y=CorrectRates,main="Rate that correct model was selected",
     xlab="Sample Size",ylab="% of Time correct model chosen",pch=16,ylim = c(0,1))
