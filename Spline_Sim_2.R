#####################################################################
#                       R Ad Hoc Functions
#####################################################################

# a general function for computing the DGOF
DGOFGenComp <- function(ICScores)
{
  DGOF <- min(ICScores)-min(ICScores[-which.min(ICScores)])
  return(DGOF)
}

# function to generate lm data
lmGen <- function(lmFit,N,NGen)
{
  # compute \hat{\sigma} for current lm fit
  tempSig <- sqrt(sum(lmFit$residuals^2)/N)
  # get the fitted values for current lm 
  tempMean <- lmFit$fitted.values
  # add tempMean to a matrix of N(0,tempSig) noise
  simDat <- tempMean + matrix(rnorm(N*NGen,mean=0,sd=tempSig),nrow=N,ncol=NGen) 
  colnames(simDat) <- paste("Draw",1:NGen,sep="")
  return(simDat)
}

#####################################################################
#             Load Packages, Rcpp Script, and Data Generation
#####################################################################

library(Rcpp)
library(RcppArmadillo)
library(splines)
# load C++ script that is used for simulation steps in ECIC
Rcpp::sourceCpp("C:\\Users\\Djizi\\OneDrive\\Desktop\\Sterner_Project\\ECIC_Code\\Spline_Sim_1_RcppCode.cpp")

# specify the interval over which the splines will be simulated
lowLim <- -10
upLim <- 10
# vector of saturated time points to plot the true spline and alternative models
manyX <- seq(from=lowLim,to=upLim,length.out=50000)
# true knots will be placed at the quintiles
trueKnots <- quantile(-10:10,probs=c(1/5,2/5,3/5,4/5))
# create a B-spline basis to plot the true spline
manyBase <- bs(manyX,knots=trueKnots,degree=3)
# set the true regression coefficients
trueCoef <- c(1,-1,1,-1,1,-1,1)
# evaluate the true spline's response values
manyY <- manyBase%*%trueCoef
# plot the true spline
plot(manyX,manyY,type="l",main=paste("Spline Plots \n
       Black=True,Red=altM1,Blue=altM2,Purple=altM3"),
     xlab="x",ylab="y")
# set up bases for alternative models
altM1Knots <- quantile(-10:10,probs=1:7/8) # knot at octiles
altM2Knots <- quantile(-10:10,probs=1:9/10) # knots at deciles
altM3Knots <- quantile(-10:10,probs=1:11/12) # knots at twelvetiles?
# create model set of knot locations for the true and alternative models
# NOTE: we fix the knots locations & degree
# but not the regression coefficients for the alternative models
M <- list("trueKnots"=trueKnots,"altM1Knots"=altM1Knots,"altM2Knots"=altM2Knots,
          "altM3Knots"=altM3Knots)
MLen <- length(M)
# fit lm using alternative bases to manyY and manyX to get an idea of how these
# models compare
lmaltM1 <- lm(manyY~bs(manyX,knots=altM1Knots,degree=3)-1)
lmaltM2 <- lm(manyY~bs(manyX,knots=altM2Knots,degree=3)-1)
lmaltM3 <- lm(manyY~bs(manyX,knots=altM3Knots,degree=3)-1)
lines(manyX,lmaltM1$fitted.values,col="red")
lines(manyX,lmaltM2$fitted.values,col="blue")
lines(manyX,lmaltM3$fitted.values,col="purple")
rm(lmaltM1)
rm(lmaltM2)
rm(lmaltM3)
rm(manyY)
rm(manyBase)
gc()
# set different sample sizes
ns <- c(15,18,25,30)
nsLen <- length(ns)
# create a list to store the observed points along the x axis for each sample size
xPoints <- list()
basisList <- list()
for(i in 1:nsLen)
{
  basisList[[i]] <- list()
}

for(i in 1:nsLen)
{
  xPoints[[i]] <- seq(from=lowLim,to=upLim,length.out=ns[i])
}

for(i in 1:nsLen)
{
  for(j in 1:MLen)
  {
    basisList[[i]][[j]] <- bs(xPoints[[i]],knots=M[[j]],degree=3)
  }
}

# set the number of draws for each sample size 
noDraws <- 500
datList <- list()
# noise added to generate data
sig=0.06
set.seed(222)
for(i in 1:nsLen)
{
  tempXPoints <- xPoints[[i]]
  tempN <- ns[i]
  tempMat <- matrix(NA,nrow=tempN,ncol=noDraws)
  tempBasisVals <- bs(tempXPoints,knots=trueKnots,degree=3)
  tempYVals <- tempBasisVals%*%trueCoef
  for(j in 1:noDraws)
  {
    draw <- tempYVals + rnorm(tempN,0,sig)
    tempMat[,j] <- draw
  }
  colnames(tempMat) <- paste("Draw",1:noDraws,sep="")
  datList[[i]] <- tempMat
}
names(datList) <- paste("Draws of n=",ns,sep="")

#####################################################################
#             Begin ECIC
#####################################################################

# ECIC step 1
# fit each of the lms for each model to the data sets
lmFitList <- list()
for(i in 1:nsLen)
{
  lmFitList[[i]] <- list()
}

for(j in 1:MLen)
{
  lmFitList[[i]][[j]] <- list()
}

for(i in 1:nsLen) # i indexes sample size
{
  tempDat <- datList[[i]]
  tempX <- xPoints[[i]]
  for(j in 1:MLen) #j indexes the model
  {
    tempBasis <- bs(tempX,knots=M[[j]],degree=3)
    lmFitList[[i]][[j]] <- apply(X=tempDat,MARGIN=2,FUN=function(x) lm(x~tempBasis-1))
  }
  names(lmFitList[[i]]) <- names(M)
  print(i)
}
names(lmFitList) <- paste("Draws of n=",ns,sep="")

# now compute the BIC under each model
BICList <- list()
for(i in 1:nsLen)
{
  BICList[[i]] <- list()
}

for(i in 1:nsLen) # i indexes the sample size
{
  tempMat <- matrix(data=NA,nrow=MLen,ncol=noDraws)
  rownames(tempMat) <- paste("BICs",names(M))
  for(j in 1:MLen) #j indexes the model
  {
    tempBICs <- sapply(lmFitList[[i]][[j]],FUN=function(x) BIC(x))
    tempMat[j,] <- tempBICs
  }
  BICList[[i]] <- tempMat
}
names(BICList) <- paste("Draws of n=",ns,sep="")

# ECIC step #2
# apply which.min(x) across each column of each matrix in BICList to store
# the best oberseved models and their corresponding element index in the model set M
# (the latter will be used for the bias correction later as in ECIC step 4a)
MbList <- list()
MbIndList <- list()
for(i in 1:nsLen)
{
  MbList[[i]] <- apply(X=BICList[[i]],MARGIN=2,FUN=function(x) names(M)[which.min(x)])
  MbIndList[[i]] <- apply(X=BICList[[i]],MARGIN=2,FUN=function(x) which.min(x))
}
names(MbList) <- paste("Draws of n=",ns,sep="")
names(MbIndList) <- paste("Draws of n=",ns,sep="")

# make barplots for the percentage of the time each model was chosen as best for each sample size
for(i in 1:nsLen)
{
  # get the frequencies of models selected as best
  tempTable <- table(MbList[[i]])
  barplot(tempTable/noDraws,main=paste("Dist. of Best Observed Models  (1 is true model)\n n=",ns[i]))
}

# list to store the observed best DGOF for each draw of each sample size
obsDGOFs <- list() 
# ECIC step #3
# compute the observed DGOFs
for(i in 1:nsLen)
{
  # compute the observed DGOFs for each draw of each sample size
  obsDGOFs[[i]] <- apply(BICList[[i]],MARGIN=2,FUN=function(x) DGOFGenComp(x))
}
names(obsDGOFs) <- paste("Draws of n=",ns,sep="")

rm(BICList)
gc()
# sample size for estimating the probability of choosing the observed best model under the assumption an
# alternative model is true
N1 <- 500
# sample size for simulating the DGOF distribution under the assumption that an alternative model is true
N2 <- 1000
# pre-specified type-1 error rate
alpha <- 0.15
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
  for(j in 1:MLen)
  {
    simDat1List[[i]][[j]] <- list()
    simDat2List[[i]][[j]] <- list()
  }
}

set.seed(19)
ptm=proc.time() #Start timing
# create matrices of random errors using \hat{\sigma} for each model
for(i in 1:nsLen) # index the sample sizes
{
  tempN <- ns[i]
  for(j in 1:MLen) # indexes the models in the model set
  {
    for(k in 1:noDraws) # index the draws
    {
      # compute \hat{\sigma} for current lm fit
      tempSig <- sqrt(sum(lmFitList[[i]][[j]][[k]]$residuals^2)/tempN) 
      # get the fitted values for current lm 
      tempMean <- lmFitList[[i]][[j]][[k]]$fitted.values
      # add tempMean to a matrix of N(0,tempSig) noise
      simDat1List[[i]][[j]][[k]] <- tempMean + matrix(rnorm(tempN*N1,mean=0,sd=tempSig),nrow=tempN,ncol=N1) 
      simDat2List[[i]][[j]][[k]] <- tempMean + matrix(rnorm(tempN*N2,mean=0,sd=tempSig),nrow=tempN,ncol=N2)
      colnames(simDat1List[[i]][[j]][[k]]) <- paste("Draw",1:N1,sep="")
      colnames(simDat2List[[i]][[j]][[k]]) <- paste("Draw",1:N2,sep="")
    }
    names(simDat1List[[i]][[j]]) <- paste("lm Fit for Obs",1:noDraws)
    names(simDat2List[[i]][[j]]) <- paste("lm Fit for Obs",1:noDraws)
  }
  names(simDat1List[[i]]) <- paste("Generated from",names(M))
  names(simDat2List[[i]]) <- paste("Generated from",names(M))
  print(i)
}
names(simDat1List) <- paste("Draws of n=",ns,sep="")
names(simDat2List) <- paste("Draws of n=",ns,sep="")
proc.time() - ptm
# create list of bases to simulate data
basisMats <- list()
for(i in 1:nsLen)
{
  basisMats[[i]] <- list()
}

for(i in 1:nsLen)
{
  for(j in 1:MLen)
  {
    basisMats[[i]][[j]] <- matrix(basisList[[i]][[j]],nrow=nrow(basisList[[i]][[j]]),
                                  ncol=ncol(basisList[[i]][[j]]))
  }
  names(basisMats[[i]]) <- paste("Basis Using",names(M))
}
names(basisMats) <- paste("Draws of n=",ns,sep="")

# Use C++ code to simulate datasets
# the RcpplmComps inputs the list of simulated draws from the MLE regression fits 
# to the observations and a list of the design matrices derived from the model set M
# matrices of BICs are returned 
set.seed(1000)
ptm=proc.time() #Start timing 
ICsSimDat1 <- RcpplmComps(simDat1List,basisMats)
proc.time() - ptm

#library(foreach)
#library(doParallel)
# put function in package and then export it into foreach?
#nCores <- parallel::detectCores()-1
# two options for parallel backends:PSOCK and FORK, FORK is faster but doesn't work on Windows :(
#myCluster <- parallel::makeCluster(
#  nCores, 
#  type = "PSOCK"
#)
#print(myCluster)
#doParallel::registerDoParallel(cl=myCluster)  # use multicore, set to the number of our cores
#check if it is registered 
#foreach::getDoParRegistered()
#how many workers are available? 
#foreach::getDoParWorkers()
#test <- foreach (i=1:4,.packages = "Rcpp",.noexport = "RcpplmCompsParallel") %dopar% 
#{
#  RcpplmCompsParallel(simDat1List[[i]],basisMats[[i]])
#}
#parallel::stopCluster(cl=myCluster)

# label the elements in ICsSimDat1
for(i in 1:nsLen)
{
  for(j in 1:MLen)
  {
    for(k in 1:noDraws)
    {
      colnames(ICsSimDat1[[i]][[j]][[k]]) <- paste("BIC Under", names(M))
    }
    names(ICsSimDat1[[i]][[j]]) <- paste("lm Fit for Obs",1:noDraws)
  }
  names(ICsSimDat1[[i]]) <- paste("Generated from",names(M))
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
      minICList[[i]][[j]][[k]] <- apply(ICsSimDat1[[i]][[j]][[k]],MARGIN=1,FUN=function(x) names(M)[which.min(x)])
    }
    names(minICList[[i]][[j]]) <- paste("lm Fit for Obs",1:noDraws)
  }
  names(minICList[[i]]) <- paste("Generated from",names(M))
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
  names(piHatList[[i]]) <- paste("lm Fits for Obs",1:noDraws)
}
names(piHatList) <- names(piHatList) <- paste("n=",ns,sep="")

for(i in 1:nsLen) # indexes the sample size
{
  for(k in 1:noDraws) # indexes the times that the kth prob is selected as best
  {
    tempMat <- matrix(data=NA,nrow=MLen,ncol=MLen)
    rownames(tempMat) <- paste("true knots=",names(M),sep="")
    colnames(tempMat) <- paste("% of time knots=",names(M)," chosen ",sep="")
    for(d in 1:MLen) # index the true knots
      for(f in 1:MLen) # index the knots being considered for modeling
        tempMat[d,f] <- sum(minICList[[i]][[d]][[k]]==names(M)[f])/N1
    piHatList[[i]][[k]] <- tempMat
  }
  print(i)
}
names(piHatList) <- paste("n=",ns,sep="")

rm(ICsSimDat1)
rm(minICList)
gc()
# list to store the quantiles from alpha/P(Mb=M)
#tauHatList <- list()
#for(i in 1:nsLen)
#{
#  tauHatList[[i]] <- list()
#  for(j in 1:MLen)
#    tauHatList[[i]][[j]] <- list()
#}

#for(i in 1:nsLen) # indexes the sample size
#{
#  for(j in 1:MLen) # indexes the true generating probability
#  {
#    for(k in 1:noDraws) # indexes the times that the kth prob is selected as best
#    {
#      tempMat <- piHatList[[i]][[j]][[k]]
#      tempMat <- alpha/tempMat
# if alpha/pi_hat >1, adjusted the value down to 1
#      tempMat[tempMat>1] <- 1
#      tauHatList[[i]][[j]][[k]] <- tempMat
#    }
#  }
#}
#names(tauHatList) <- paste("n=",ns,sep="")

# Use C++ code to simulate datasets
# the RcppSimDGOFs inputs the list of simulated draws from the MLE regression fits 
# to the observations and a list of the design matrices derived from the model set M
# DGOF distributions are returned
ptm=proc.time() #Start timing 
DGOFList <- RcppSimDGOFs(simDat2List,basisMats)
proc.time() - ptm
# label the elements in DGOFList
for(i in 1:nsLen)
{
  for(j in 1:MLen)
  {
    for(k in 1:noDraws)
    {
      colnames(DGOFList[[i]][[j]][[k]]) <- paste("DGOF Under", names(M)," Observed Best")
    }
    names(DGOFList[[i]][[j]]) <- paste("lm Fit for Obs",1:noDraws)
  }
  names(DGOFList[[i]]) <- paste("Generated from",names(M))
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
    tempMbInd <- which(names(M)==MbList[[i]][k])
    # identify the observed DGOF
    tempObsDGOF <- obsDGOFs[[i]][k]
    # identify the alternative models
    altModels <- names(M)[-tempMbInd]
    tempDGOFQuantiles <- rep(NA,MLen-1)
    for(h in 1:(MLen-1)) # assume h is index for the true model
    {
      # now just retrieve the quantile and DGOF distribution
      # to compare it to the observed DGOF
      curAltModel <- altModels[h]
      curAltModelInd <- which(names(M)==curAltModel)
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
names(assessList) <- names(DGOFList) <- paste("Draws of n=",ns,sep="")

DecisionRates <- sapply(X=assessList,FUN=function(x) sum(as.numeric(x[2,]))/noDraws)
plot(x=ns,y=DecisionRates,main="Proportion of runs a model was selected",
     xlab="Sample Size",ylab="% of a model was selected",pch=16,ylim = c(0,1))
# take a look at type 1 error rate
# subset assess list by only when a decision was made i.e. second row = 1
decAssessList <- lapply(X=assessList,FUN=function(x) x[,x[2,]=="1"])
# compute type 1 error rate
t1ErrorRates <- sapply(X=decAssessList,FUN=function(x) sum(x[1,]!="trueKnots")/noDraws)
plot(x=ns,y=t1ErrorRates,main=c("Type 1 Error Rates by Sample Size at alpha=",alpha),
     xlab="Sample Size",ylab="% of Time Wrong Best Model Selected",pch=16,ylim = c(0,2*alpha))
CorrectRates <- sapply(X=decAssessList,FUN=function(x) sum(x[1,]=="trueKnots")/noDraws)
plot(x=ns,y=CorrectRates,main="Rate that correct model was selected",
     xlab="Sample Size",ylab="% of Time correct model chosen",pch=16,ylim = c(0,1))


