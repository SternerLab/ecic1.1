# compute the log likelihood for a Bernoulli distribution 
LLBern <- function(dat,p)
{
  LL <- log(p)*sum(dat) + log(1-p)*sum(1-dat)
  return(LL)
}

# compute the log likelihood for a normal distribution with fixed sigma
# LL=-n/2*log(2\pi)-n/2*log(\sigma^2)-1/(2\sigma^2)\sum(x_j-\mu)^2
LLNorm <- function(dat,pars)
{
  mu=pars[1]
  sig2=pars[2]^2
  n=length(dat)
  # since sigma is fixed, we only need the third term in the LL
  LL <- -n/2*log(2*pi)-n/2*log(sig2)-1/(2*sig2)*sum((dat-mu)^2)
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

# generate a list of random Bernoulli draws where each element holds
# a matrix with noDraws draws of size ns[i] 
generateData <- function(tempN,noDraws,p,dataType)
{
  if(dataType=="Bernoulli")
    tempDraws <- rbinom(n=tempN*noDraws,size=1,prob=p)
  else
    tempDraws <- rnorm(n=tempN*noDraws,mean=p[1],sd=p[2])
  # convert vector of draws into a matrix of multiple draws
  tempDraws <- matrix(data=tempDraws,nrow=tempN,ncol=noDraws)
  # store matrix as the draws that are each of sample size tempN
  colnames(tempDraws) <- paste("Draw",1:noDraws,sep="")
  return(tempDraws)
}

# compute the IC under each model in the model set for each sample size
ICComputations <- function(datMat,M,MLen,noDraws,ICType)
{
  # MLenXnoDraws matrix holding the IC value for the 3 models for each draw
  tempMat <- matrix(NA,nrow=MLen,ncol=noDraws)
  if(ICType=="BernoulliNegLL")
  {
    for(j in 1:MLen)
    {
      tempICs <- -1*apply(X=datMat,MARGIN=2,FUN=function(x) LLBern(x,p=M[j]))
      tempMat[j,] <- tempICs
    }
    rownames(tempMat) <- paste("-LL Under p=",M,sep="")
    colnames(tempMat) <- paste("Draw",1:noDraws,sep="")
  }
  else
  {
    for(j in 1:MLen)
    {
      tempICs <- -1*apply(X=datMat,MARGIN=2,FUN=function(x) LLNorm(x,p=M[[j]]))
      tempMat[j,] <- tempICs
    }
    rownames(tempMat) <- paste("-LL Under N(0",M,sep="")
    colnames(tempMat) <- paste("Draw",1:noDraws,sep="")
  }
  return(tempMat)
}

# apply which.min(x) across each column of each matrix in the ICComps list 
# to determine the observed best models
# list of length nsLen of 3XnoDraws vectors in each element that holds the lowest 
# IC score for each draw of each sample size
MbComputations <- function(ICCompsMat,MNames)
{
  Mbs <- apply(X=ICCompsMat,MARGIN=2,FUN=function(x) MNames[which.min(x)])
  return(Mbs)
}

obsDGOFsComputations <- function(ICComps,nsLen)
{
  obsDGOFs <- list()
  for(i in 1:nsLen)
  {
    # compute the observed DGOFs for each draw of each sample size
    obsDGOFs[[i]] <- apply(X=ICComps[[i]],MARGIN=2,FUN=function(x) DGOFGenComp(x))
  }
  names(obsDGOFs) <- paste("Observed DGOFs for ","Draws of n=",ns,sep="")
  return(obsDGOFs)
}

# compute matrices that hold P_i(g(F)=M_b)
piHatMatComputations <- function(minICListFixedN,MLen,N1,MNames)
{
  tempMat <- matrix(data=NA,nrow=MLen,ncol=MLen)
  rownames(tempMat) <- paste("true p=",M,sep="")
  colnames(tempMat) <- paste("% of time p=",M," observed best ",sep="")
  for(j in 1:MLen) # indexes the true generating probability
  {
    for(k in 1:MLen) # indexes the times that the kth prob is selected as best
    {
      tempMat[j,k] <- sum(minICList[[i]][[j]]==MNames[k])/N1
    }
  }
  return(tempMat)
}

DGOFSimComputations <- function(ICsSimDat2,nsLen,mLen)
{
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
  return(DGOFList)
}
# go through each of the alternative models, assume they are true, and compute quantiles
# i indexes sample size, j indexes models in the alternative set, and k indexes models in the full set
# list to store accept or reject decisions by sample size
ECICDecisions <- function(MbList,obsDGOFs,piHatList,DGOFList,alpha,M,nsLen,MLen,noDraws)
{
  thresholds <- list()
  aOrRList <- list()
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
        probThisFalseMbInd <- piHatList[[i]][curAltModelInd,tempMbInd]
        if(probThisFalseMbInd==0) 
          tempTau <- 1
        else
        {
          totalProbFalseMb <- sum(piHatList[[i]][curAltModelInd,-curAltModelInd])
          tempTau <- alpha*probThisFalseMbInd/totalProbFalseMb
          if(tempTau>1) 
            tempTau <- 1
        }
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
    aOrRList[[i]] <- tempAorRVec
  }
  returnList <- list("thresholds"=thresholds,"aorRList"=aOrRList)
}
