#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
#define pi 3.14159265358979
// [[Rcpp::depends(RcppArmadillo)]]

// auxiliary function
// [[Rcpp::export]]
arma::mat DGOFCompforMatrix(arma::mat ICScores)
{
    int ncols = ICScores.n_cols;
    int nrows = ICScores.n_rows;
    arma::mat returnMat(nrows,ncols);
    for(int i=0;i<ncols;++i)
    {
        arma::mat tempSubMat = ICScores;
        tempSubMat.shed_col(i);
        arma::vec tempFillVec(nrows);
        arma::vec tempMbVec = ICScores.col(i);
        for(int j=0;j<nrows;++j)
        {
            arma::vec tempRow = vectorise(tempSubMat.row(j));
            double tempDGOF = tempMbVec[j] - tempRow.min();
            tempFillVec[j] = tempDGOF;
        }
        returnMat.col(i) = tempFillVec;
    }
    return(returnMat);
}


// compute the BIC for a normally distributed random variable
// [[Rcpp::export]]
double normBICComp(arma::vec ys,arma::vec meanEst, double varEst,double p)
{
    double N = ys.n_elem;
    arma::vec sqrdDiffs = pow(ys-meanEst,2);
    double LL = -N/2.0*log(2.0*pi)-N/2.0*log(varEst)-accu(sqrdDiffs)/(2.0*varEst);
    double BIC = -2.0*LL+p*log(N);
    return(BIC);
}

// compute the BIC for draws in a data frame. Returns a vector of BICs
// [[Rcpp::export]]
arma::vec normBICComp2(arma::mat tempYMat,arma::mat tempHatMat,int noObsPts,int p,int noSimDraws,arma::vec fillVec)
{
    for(int d=0;d<noSimDraws;++d) // indexes data generated from a specific regression fit
    {
        arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
        arma::vec tempMean = tempHatMat*tempYs; // estimate the regression coefficients
        double tempVar = (accu(pow(tempYs-tempMean,2)))/noObsPts; // estimate the population sd under MLE
        // Compute the BIC
        arma::vec sqrdDiffs = pow(tempYs-tempMean,2);
        double LL = -noObsPts/2.0*log(2.0*pi)-noObsPts/2.0*log(tempVar)-accu(sqrdDiffs)/(2.0*tempVar);
        fillVec(d) = -2.0*LL+p*log(noObsPts);
    }
    return(fillVec);
}

// compute the BIC for draws in a data frame. Returns a vector of BICs
// includes the MLE adjustment
// [[Rcpp::export]]
arma::vec normBICComp3(List tempListX1,arma::mat tempYMat,arma::mat tempX, arma::mat invR2Q2,arma::vec fillVec,
                       int tempMbInd,int p,int noObsPts,int noSimDraws,int mLen)
{
    arma::mat tempHatMat=tempX*invR2Q2;
    for(int d=0;d<noSimDraws;++d) // indexes data generated from a specific regression fit
    {
        arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
        arma::vec tempBetaHat=invR2Q2*tempYs;
        arma::vec tempMean = tempHatMat*tempYs; // estimate the regression coefficients
         arma::vec sqrdDiffsNotCor = pow(tempYs-tempMean,2);
        double tempVar = accu(sqrdDiffsNotCor)/noObsPts; // estimate the population sd under MLE
        // perform bias correction on both tempBetaHat and tempVar
        int itCount = 0; // counter for the number of iterations to try sampling new data sets as in step 4a2 in ECIC
        int MbObsCount = 0; // counter for the number of times Mb is observed best using the new datasets
        double varSum = 0; // accumulator used to adjust the MLE variance later
        arma::vec coeffSum = zeros(p); // accumulator used to adjust the MLE coefficients later
        int stopIts = 1500; // set a threshold to give up on trying to generate N1 new datasets for step 4a2 in ECIC
        int sufficientSamps = 30; // set a sufficient number of times we need to observe the given model as best
        arma::mat mvnCov = eye(noObsPts,noObsPts);
        mvnCov = mvnCov*tempVar;
        while(itCount<stopIts && MbObsCount<sufficientSamps) // stop whenever the "giving up threshold" is reached or N1 samples where Mb is observed best occur
        {
            arma::vec newDraw = mvnrnd(tempMean,mvnCov); // draw a new data set using the MLEs
            arma::vec modelBICs = vec(mLen); // create a vector to hold the BICs for each model under the newDraw
            for(int j=0;j<mLen;++j) // go through each model in the model set, and see which one has the lowest BIC using the previously generated data
            {
                arma::mat tempX2 = tempListX1[j]; // retrieve the design matrix for the current model in the model set
                int ncols = tempX2.n_cols; // number of columns in design matrix X i.e. the number of regression coefficients to estimate
                double p2 = ncols+1; // number of total parameters to be estimated (regression coefficients + 1 for the variance)
                // compute a QR decomposition of tempX for more stable computations
                arma::mat Q,R;
                // perform the decomposition
                arma::qr(Q,R,tempX);
                // subset relevant parts of the Q and R matrices
                arma::mat R2 = R.rows(0,ncols-1);
                arma::mat Q2 = Q.cols(0,ncols-1);
                // compute the beta hat
                arma::vec tempBeta2 = inv(R2)*Q2.t()*newDraw;
                // compute the expectation function
                arma::vec tempMean2 = tempX2*tempBeta2;
                arma::vec sqrdDiffs2 = pow(newDraw-tempMean2,2);
                double tempVar2 = accu(sqrdDiffs2)/noObsPts; // MLE variance
                double LL = -noObsPts/2.0*log(2.0*pi)-noObsPts/2.0*log(tempVar2)-accu(sqrdDiffs2)/(2.0*tempVar2);
                modelBICs(j) = -2.0*LL+p2*log(noObsPts);
            }
            // find model with the minimum BIC
            int minBICModel = index_min(modelBICs);
            if(minBICModel==tempMbInd) // if the model with minimum BIC is Mb then add to accumulator to be averaged later
            {
                MbObsCount++;
                // compute the MLE of the current model under the current simulated data set
                arma::vec newCoeffs = invR2Q2*newDraw;
                arma::vec newMean = tempX*newCoeffs;
                double newVar = (accu(pow(newDraw-newMean,2)))/noObsPts; // estimate the population var under MLE
                varSum+=newVar;
                coeffSum+=newCoeffs;
            }
            itCount++;
        }
    // If sufficientSamps where M_b is observed as best are found then return the bias corrected parameters. Otherwise return the unadjusted MLEs
        if(MbObsCount==sufficientSamps)
        {
            arma::vec biasBetaHat = coeffSum/sufficientSamps - tempBetaHat;
            double biasVar = varSum/sufficientSamps-tempVar;
            arma::vec biasCorrBetaHat = tempBetaHat-biasBetaHat;
            double biasCorrVar = tempVar-biasVar;
            arma::vec biasCorrMean = tempX*biasCorrBetaHat;
            // Compute the BIC
            arma::vec sqrdDiffsCor = pow(tempYs-biasCorrMean,2);
            double LL = -noObsPts/2.0*log(2.0*pi)-noObsPts/2.0*log(biasCorrVar)-accu(sqrdDiffsCor)/(2.0*biasCorrVar);
            fillVec(d) = -2.0*LL+p*log(noObsPts);
        }
        else
        {
            // Compute the BIC
            double LL = -noObsPts/2.0*log(2.0*pi)-noObsPts/2.0*log(tempVar)-accu(sqrdDiffsNotCor)/(2.0*tempVar);
            fillVec(d) = -2.0*LL+p*log(noObsPts);
        }
    }
    return(fillVec);
}

// [[Rcpp::export]]
List lmCompsRcpp(List yList,List basisList)
{
    int nLen = basisList.length(); // number of sample sizes
    List returnList(nLen); // the return list that indexes the sample sizes
    for(int i=0;i<nLen;++i) // indexes sample size
    {
        List tempListX1 = basisList[i]; // step into design matrices X of a particular sample size
        int mLen = tempListX1.length(); // number of models
        List returnListSub1(mLen); // part of return list that indexes the generating model
        List tempListY1 = yList[i]; // step into observations of a particular sample size
        for(int j=0;j<mLen;++j) // indexes true generating model
        {
            List tempListY2 = tempListY1[j]; // step into observations generated from a particular model
            int noDraws = tempListY2.length(); // number of draws
            List returnListSub2(noDraws); // part return list that indexes the exact fits for each model
            for(int k=0;k<noDraws;++k) // indexes specific regression fit
            {
                arma::mat tempYMat = tempListY2[k]; // current simulated draws for sample size i, true model j, and model fit k
                int noSimDraws = tempYMat.n_cols; // (N1 in ECIC paper notation)
                // create matrix to store BICs for data generated from a specific model fit
                arma::mat fillMat(noSimDraws,mLen);
                int noObsPts = tempYMat.n_rows; // number of observed points for each simulated draw
                for(int m=0;m<mLen;++m) // indexes regression model to fit for BIC
                {
                    arma::mat tempX = tempListX1[m]; // design matrix X for current sample size and model
                    int ncols = tempX.n_cols; // number of columns in design matrix X i.e. the number of regression coefficients to estimate
                    double p = ncols+1; // number of total parameters to be estimated (regression coefficients + 1 for the population sd)
                    // compute a QR decomposition of X for more stable hat matrix computation
                    arma::mat Q,R;
                    // perform the decomposition
                    arma::qr(Q,R,tempX);
                    // subset relevant parts of the Q and R matrices
                    arma::mat R2 = R.rows(0,ncols-1);
                    arma::mat Q2 = Q.cols(0,ncols-1);
                    // compute the hat matrix
                    arma::mat tempHatMat = tempX*inv(R2)*Q2.t();
                    // vector to fill the BIC computation for each simulated draw generated from model j but assuming here that model m is true
                    arma::vec fillVec = vec(noSimDraws);
                    fillVec=normBICComp2(tempYMat,tempHatMat,noObsPts,p,noSimDraws,fillVec);
                    fillMat.col(m) = fillVec; // column m of fillMat will hold the BICs for each draw computed under model m
                }
                returnListSub2[k] = fillMat; // store all BIC computations for regression fit k
            }
        returnListSub1[j] = returnListSub2; // store the list of BIC computations for all noDraws matrices truly generated from model j
        }
    returnList[i] = returnListSub1; // store the list of BIC computations for all noDraws matrices truly generated from model j with sample sizes i
    Rcout << "Its. for Sample Size Index " << i+1 << " completed" << endl;
    }
    return(returnList);
}


// [[Rcpp::export]]
List RcpplmCompswMLEAdjust(List yList,List basisList,List MbIndList)
{
    int nLen = basisList.length(); // number of sample sizes
    List returnList(nLen); // the return list that indexes the sample sizes
    for(int i=0;i<nLen;i++) // indexes sample size
    {
        arma::vec tempMbs = MbIndList[i]; // step into the best observed models of a particular sample size
        List tempListX1 = basisList[i]; // step into design matrices X of a particular sample size
        int mLen = tempListX1.length(); // number of models
        List returnListSub1(mLen); // part of return list that indexes the generating model
        List tempListY1 = yList[i]; // step into observations of a particular sample size
        for(int j=0;j<mLen;j++) // indexes true generating model
        {
            List tempListY2 = tempListY1[j]; // step into observations generated from a particular model
            int noDraws = tempListY2.length(); // number of draws
            List returnListSub2(noDraws); // part return list that indexes the exact fits for each model
            for(int k=0;k<noDraws;k++) // indexes specific regression fit
            {
                int tempMbIndCVersion = tempMbs[k]-1; // retrieve the index for the best observed model for the draw (subtract 1 to convert R index into C++ index)
                arma::mat tempYMat = tempListY2[k]; // current simulated draws for sample size i, true model j, and model fit k
                int noSimDraws = tempYMat.n_cols; // (N1 in ECIC paper notation)
                // create matrix to store BICs for data generated from a specific model fit
                arma::mat fillMat(noSimDraws,mLen);
                int noObsPts = tempYMat.n_rows; // number of observed points for each simulated draw
                for(int m=0;m<mLen;m++) // indexes regression model to fit for BIC
                {
                    arma::mat tempX = tempListX1[m]; // design matrix X for current sample size and model
                    int ncols = tempX.n_cols; // number of columns in design matrix X i.e. the number of regression coefficients to estimate
                    double p = ncols+1; // number of total parameters to be estimated (regression coefficients + 1 for the population sd)
                    // compute a QR decomposition of X for more stable hat matrix computation
                    arma::mat Q,R;
                    // perform the decomposition
                    arma::qr(Q,R,tempX);
                    // subset relevant parts of the Q and R matrices
                    arma::mat R2 = R.rows(0,ncols-1);
                    arma::mat Q2 = Q.cols(0,ncols-1);
                    arma::mat invR2Q2 = inv(R2)*Q2.t();
                    // vector to fill the BIC computation for each simulated draw generated from model j but assuming here that model m is true
                    arma::vec fillVec = vec(noSimDraws);
                    fillVec = normBICComp3(tempListX1,tempYMat,tempX,invR2Q2,fillVec,tempMbIndCVersion,p,noObsPts,noSimDraws,mLen); // compute the BIC for this fit
                    fillMat.col(m) = fillVec; // column m of fillMat will hold the BICs for each draw computed under model m
                }
                returnListSub2[k] = fillMat; // store all BIC computations for regression fit k;
            }
        returnListSub1[j] = returnListSub2; // store the list of BIC computations for all noDraws matrices truly generated from model j
        Rcout << "Model " << j+1 << " sample size " << i+1 << " completed" << endl;
        }
    returnList[i] = returnListSub1; // store the list of BIC computations for all noDraws matrices truly generated from model j with sample sizes i
    Rcout << "Its. for Sample Size Index " << i+1 << " completed" << endl;
    }
    return(returnList);
}

// pre-computes the DGOF distributions (code is more or less the same as RcpplmComps except for the arma::mat tempDGOFs = DGOFCompforMatrix(fillMat) line)
// [[Rcpp::export]]
List lmDGOFsRcpp(List yList,List basisList)
{
    int nLen = basisList.length(); // number of sample sizes
    List returnList(nLen); // return list that indexes the sample sizes
    for(int i=0;i<nLen;++i) // indexes sample size
    {
        List tempListX1 = basisList[i]; // step into design matrices X of a particular sample size
        int mLen = tempListX1.length(); // number of models
        List returnListSub1(mLen); // part of return list that indexes the generating model
        List tempListY1 = yList[i]; // step into observations of a particular sample size
        for(int j=0;j<mLen;++j) // indexes true generating model
        {
            List tempListY2 = tempListY1[j]; // step into observations generated from a particular model
            int noDraws = tempListY2.length(); // number of draws
            List returnListSub2(noDraws); // part of return list that indexes the exact fits for each model
            for(int k=0;k<noDraws;k++) // indexes specific regression fit
            {
                arma::mat tempYMat = tempListY2[k]; // current simulation distribution of values
                int noSimDraws = tempYMat.n_cols; // (N2 in ECIC paper notation)
                // create matrix to store BICs for data generated from a specific model fit
                arma::mat fillMat(noSimDraws,mLen);
                int noObsPts = tempYMat.n_rows; // number of observed points for each simulated draw
                for(int m=0;m<mLen;++m) // indexes regression model to fit for BIC
                {
                    arma::mat tempX = tempListX1[m]; // X for current sample size and model
                    int ncols = tempX.n_cols; // number of columns in design matrix X i.e. the number of regression coefficients to estimate
                    double p = ncols+1; // number of total parameters to be estimated (regression coefficients + 1 for the population sd)
                    // compute a QR decomposition of X for more stable hat matrix computation
                    arma::mat Q,R;
                    // perform the decomposition
                    arma::qr(Q,R,tempX);
                    // subset relevant parts of the Q and R matrices
                    arma::mat R2 = R.rows(0,ncols-1);
                    arma::mat Q2 = Q.cols(0,ncols-1);
                    // compute the hat matrix
                    arma::mat tempHatMat = tempX*inv(R2)*Q2.t();
                    // vector to fill the BIC computation for each simulated draw generated from model j but assuming here that model m is true
                    arma::vec fillVec = vec(noSimDraws);
                    fillVec=normBICComp2(tempYMat,tempHatMat,noObsPts,p,noSimDraws,fillVec);
                    fillMat.col(m) = fillVec;
                }
                arma::mat tempDGOFs = DGOFCompforMatrix(fillMat);
                returnListSub2[k] = tempDGOFs;
            }
        returnListSub1[j] = returnListSub2;
        }
    returnList[i] = returnListSub1;
    Rcout << "Its. for Sample Size Index " << i+1 << " completed" << endl;
    }
    return(returnList);
}

