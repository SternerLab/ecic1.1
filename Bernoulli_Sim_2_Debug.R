#####################################################################
#                       Description
#####################################################################

# This script assesses error rates
# for a model set of three Bernoulli distributions with parameters
# 0.48,0.65, and 0.75
# We simulate a scenario where the true probability distribution is a Bernoulli 
# distribution with p=0.65
# The purpose of this simulation study is debug the code in 
# Bernoulli_Sim_2.R

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
# store the true generating probability
trueP <- 0.65
# set the sample size for each draw
n <- 20
# set the number of draws 
noDraws <- 2000
set.seed(222)
# simulate data 
# draw from a Bernoulli dist with saved n & p
tempDraws <- rbinom(n=n*noDraws,size=1,prob=trueP)
# convert vector of draws into a matrix of multiple draws
dat <- matrix(data=tempDraws,nrow=n,ncol=noDraws)
colnames(dat) <- paste("Draw",1:noDraws,sep="")

#####################################################################
#                      Run Study for Scenario 1
#####################################################################

# ECIC step #1
# compute the IC under each model in the model set for each sample size
ICComps <- matrix(NA,nrow=MLen,ncol=noDraws)
for(j in 1:MLen)
{
  tempICs <- -1*apply(X=dat,MARGIN=2,FUN=function(x) LLBern(x,p=M[j]))
  ICComps[j,] <- tempICs
}
# columns are draws and rows are the IC evaluated under different parameter values
rownames(ICComps) <- paste("IC Under p=",M,sep="")
colnames(ICComps) <- paste("Draw",1:noDraws,sep="")

# ECIC step #2
# apply which.min(x) across each column in ICComps 
# to determine the observed best models
# note 1 corresponds to p=0.5, 2 to p=0.65, and 3 to 0.75
MbVec <- apply(X=ICComps,MARGIN=2,FUN=function(x) M[which.min(x)])

# make barplot for the percentage of the time each model was chosen as best
tempTable <- table(MbVec)
barplot(tempTable/noDraws,main=paste("Dist. of Best Observed Models  (p=0.65 is true model)\n n=",n))

# ECIC step #3
# compute the observed DGOFs
obsDGOFs <- apply(X=ICComps,MARGIN=2,FUN=function(x) DGOFGenComp(x))
hist(obsDGOFs,main="Distribution of the Observed DGOFs")
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
set.seed(19)
# simulate data for each probability in the model set
for(j in 1:MLen) #  j indexes the assumed true probability
{
  simDat1List[[j]] <- rbinom(n=n*N1,size=1,prob=M[j])
  simDat1List[[j]] <- matrix(data=simDat1List[[j]],nrow=n,ncol=N1)
  colnames(simDat1List[[j]]) <- paste("Draw",1:N1,sep="")
  simDat2List[[j]] <- rbinom(n=n*N2,size=1,prob=M[j])
  simDat2List[[j]] <- matrix(data=simDat2List[[j]],nrow=n,ncol=N2)
  colnames(simDat2List[[j]]) <- paste("Draw",1:N2,sep="")
}
names(simDat1List) <- paste("True p=",M,sep="")
names(simDat2List) <- paste("True p=",M,sep="")

# each element in these lists will hold matrices of the -LL's for all draws 
# under all the probabilities in the model set 
ICsSimDat1 <- list()
ICsSimDat2 <- list()
for(j in 1:MLen) # j indexes the assumed true probability 
{
  # MLenXN1 matrix to store -LL's 
  tempMat1 <- matrix(NA,nrow=MLen,ncol=N1)
  colnames(tempMat1) <- paste("Draw",1:N1,sep="")
  # MLenXN2 matrix to store -LL's 
  tempMat2 <- matrix(NA,nrow=MLen,ncol=N2)
  colnames(tempMat2) <- paste("Draw",1:N2,sep="")
  tempDraws1 <- simDat1List[[j]] 
  tempDraws2 <- simDat2List[[j]] 
  for(k in 1:MLen) # compute -LL's under different parameters in the model set
  {
    tempICs1 <- -1*apply(X=tempDraws1,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
    tempICs2 <- -1*apply(X=tempDraws2,MARGIN=2,FUN=function(x) LLBern(x,p=M[k]))
    tempMat1[k,] <- tempICs1
    tempMat2[k,] <- tempICs2
  }
  ICsSimDat1[[j]] <- tempMat1
  ICsSimDat2[[j]] <- tempMat2
  rownames(ICsSimDat1[[j]]) <- paste("-LL Under p=",M,sep="")
  rownames(ICsSimDat2[[j]]) <- paste("-LL Under p=",M,sep="")
}
names(ICsSimDat1) <- paste("True p=",M,sep="")
names(ICsSimDat2) <- paste("True p=",M,sep="")

# determine the models with the minimum IC for each set of draws
minICList <- list()
for(j in 1:MLen) # j indexes the assumed true parameter
{
  minICList[[j]] <- apply(ICsSimDat1[[j]],MARGIN=2,FUN=function(x) M[which.min(x)])
}
names(minICList) <- paste("True p=",M,sep="")

# create a matrix that holds P_i(g(F)=M_b)
piHatMat <- matrix(data=NA,nrow=MLen,ncol=MLen)
rownames(piHatMat) <- paste("true p=",M,sep="")
colnames(piHatMat) <- paste("% of time p=",M," observed best ",sep="")
for(j in 1:MLen) # indexes the true generating probability
{
  for(k in 1:MLen) # indexes the times that the kth prob is selected as best
  {
    piHatMat[j,k] <- sum(minICList[[j]]==M[k])/N1
  }
}

# list to store the quantiles from alpha/P(Mb=M)
tauHatMat <- alpha/piHatMat
tauHatMat <- alpha/piHatMat
# if alpha/pi_hat >1, adjust the value down to 1
tauHatMat[tauHatMat>1] <- 1
colnames(tauHatMat) <- paste("tauth percentile using (p=",M," observed best) ",sep="")

# compute the DGOF distributions as in step 4c) of ECIC
DGOFList <- list()
for(j in 1:MLen) # index the true generating probability 
{
  tempDat <- ICsSimDat2[[j]]
  # matrix to store the DGOF distributions under the assumption that model k is
  # the observed best index
  tempMat <- matrix(NA,nrow=MLen,ncol=N2)
  colnames(tempMat) <- paste("Draw",1:N2,sep="")
  for(k in 1:MLen) #index the model assumed to be best
  {
    tempDGOFs <- apply(tempDat,MARGIN=2,FUN=function(x) DGOFSimComp(x,MbInd=k))
    tempMat[k,] <- tempDGOFs # store row as DGOF distribution
  }
  DGOFList[[j]] <- tempMat
  rownames(DGOFList[[j]]) <- paste("DGOFs,M_b=",M,sep="")
}
names(DGOFList) <- paste("True p=",M,sep="")
# print the number of unique values for each DGOF distribution
length(unique(DGOFList[[1]][1,]))
length(unique(DGOFList[[1]][2,]))
length(unique(DGOFList[[1]][3,]))
length(unique(DGOFList[[2]][1,]))
length(unique(DGOFList[[2]][2,]))
length(unique(DGOFList[[2]][3,]))
length(unique(DGOFList[[3]][1,]))
length(unique(DGOFList[[3]][2,]))
length(unique(DGOFList[[3]][3,]))

# plot DGOF distributions and quantiles
for(i in 1:MLen)
{
  trueGenModTemp <- M[i] #fix the true generating probability of a DGOF dist.
  for(j in 1:MLen)
  {
    bestObsModTemp <- M[j] #fix the observed best model
    tempTauRound <- round(tauHatMat[i,j],3)
    tempQuant <- quantile(DGOFList[[i]][j,],tauHatMat[i,j])
    hist(DGOFList[[i]][j,],main=paste("DGOF: True p=",trueGenModTemp,",Obs. Best p=",bestObsModTemp,"\n tau (red line)=",tempTauRound),xlab="Deltaf")
    abline(v=tempQuant,col="red")
  }
}

# ECIC steps #4 and #5
# go through each of the alternative models, assume they are true, and compute quantiles
# i indexes sample size, j indexes models in the alternative set, and k indexes models in the full set
# matrix to store the decision thresholds
thresholds <- matrix(NA,nrow=noDraws,ncol=MLen-1)
# vector to store accept or reject decisions
aOrRVec <- rep(NA,noDraws)
for(j in 1:noDraws) # go through each observation and perform ECIC
{
  
  # get the element index of the best observed model 
  tempMbInd <- which(M==tempMb)
  # identify the observed DGOF
  tempObsDGOF <- obsDGOFs[j]
  # identify the alternative models
  altModels <- M[-tempMbInd]
  tempDGOFQuantiles <- rep(NA,MLen-1)
  for(k in 1:(MLen-1)) # assume k is index for the true model
  {
    # now just retrieve the quantile and DGOF distribution
    # to compare it to the observed DGOF
    curAltModel <- altModels[k]
    curAltModelInd <- which(M==curAltModel)
    tempTau <- tauHatMat[curAltModelInd,tempMbInd]
    tempDGOFs <- DGOFList[[curAltModelInd]][tempMbInd,]
    tempDGOFQuantiles[k] <- quantile(tempDGOFs,probs=tempTau)
  }
  thresholds[j,] <- tempDGOFQuantiles
    # take the alternative model quantile estimates
  tempFinQuantile <- min(tempDGOFQuantiles)
  # ECIC step #6
  # store 1 if observed model is chosen as best and 0 otherwise
  aOrRVec[j] <- ifelse(test=tempObsDGOF<tempFinQuantile,yes=1,no=0)
}

# plot the percentage of the time no decision is made by sample size
fracNoDec <- length(aOrRVec[aOrRVec==0])/noDraws
plot(x=n,y=fracNoDec,main="Percentage of Time No Decision is Made",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

assessMat <- rbind(MbVec,aOrRVec)
# Look at percentage of times the true model was observed as best but no decision was made (false negative)
# subset assessList to only contain observations where with no decision i.e. 0
noDecAssess <- assessMat[,assessMat[2,]==0]
rightObsUndec <- sum(noDecAssess[1,]==2)/noDraws

# note that once a decision is always made (which occurs asymptotically)
# there will not be a point corresponding to the sample size on this following plot
plot(x=n,y=rightObsUndec,main="Percentage of Time True Model Observed When No Decision Made",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)
# take a look at type 1 error rate
# subset assessMat by only when a decision was made i.e. second row = 1
decAssess <- assessMat[,assessMat[2,]==1]
# compute type 1 error rate
t1ErrorRate <- sum(decAssess[1,]!=2)/noDraws
plot(x=n,y=t1ErrorRate,main="Type 1 Error Rates",
     xlab="Sample Size",ylab="% of Time No Decision Made",pch=16)

# bind together the threshold quantiles, observed DGOFs, and observed best models
checkMat <- as.data.frame(cbind(thresholds,obsDGOFs,MbVec,aOrRVec))
checkMatT1 <- checkMat[checkMat$aOrRVec==1 & checkMat$MbVec!=2,]
