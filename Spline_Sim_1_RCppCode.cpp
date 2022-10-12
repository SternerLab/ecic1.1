#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// just a check to see if the QR decomposition is working as it should
// [[Rcpp::export]]
List qrtest(arma::mat tempX)
{
    int ncols = tempX.n_cols;
    arma::mat Q;
    arma::mat R;
    arma::qr(Q,R,tempX);
    arma::mat R2 = R.rows(0,ncols-1);
    arma::mat Q2 = Q.cols(0,ncols-1);
    arma::mat test = inv(tempX.t()*tempX)*tempX.t();
    arma::mat test2 = inv(R2)*Q2.t();
    return(List::create(Named("Orig") = test, _["QR Decomp"] = test2));
}

// auxiliary function
// [[Rcpp::export]]
arma::mat DGOFCompforMatrix(arma::mat ICScores)
{
    int ncols = ICScores.n_cols;
    int nrows = ICScores.n_rows;
    arma::mat returnMat(nrows,ncols);
    for(int i=0;i<ncols;i++)
    {
        arma::mat tempSubMat = ICScores;
        tempSubMat.shed_col(i);
        arma::vec tempFillVec(nrows);
        arma::vec tempMbVec = ICScores.col(i);
        for(int j=0;j<nrows;j++)
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
    double pi=3.14159265358979;
    double N = ys.n_elem;
    arma::vec sqrdDiffs = pow(ys-meanEst,2);
    double LL = -N/2.0*log(2.0*pi)-N/2.0*log(varEst)-accu(sqrdDiffs)/(2.0*varEst);
    double BIC = -2*LL+p*log(N);
    return(BIC);
}

// [[Rcpp::export]]
List RcpplmComps(List yList,List basisList)
{
    int nLen = basisList.length(); // number of sample sizes
    List returnList(nLen); // the return list that indexes the sample sizes
    for(int i=0;i<nLen;i++) // indexes sample size
    {
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
                    // compute the hat matrix
                    arma::mat tempHatMat = tempX*inv(R2)*Q2.t();
                    // vector to fill the BIC computation for each simulated draw generated from model j but assuming here that model m is true
                    arma::vec fillVec = vec(noSimDraws);
                    for(int d=0;d<noSimDraws;d++) // indexes data generated from a specific regression fit
                    {
                        arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
                        arma::vec tempMean = tempHatMat*tempYs; // estimate the regression coefficients
                        double tempVar = (accu(pow(tempYs-tempMean,2)))/noObsPts; // estimate the population sd under MLE
                        fillVec[d] = normBICComp(tempYs,tempMean,tempVar,p); // compute the BIC for this fit
                    }
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


// pre-computes the DGOF distributions (code is more or less the same as RcpplmComps except for the arma::mat tempDGOFs = DGOFCompforMatrix(fillMat) line)
// [[Rcpp::export]]
List RcppSimDGOFs(List yList,List basisList)
{
    int nLen = basisList.length(); // number of sample sizes
    List returnList(nLen); // return list that indexes the sample sizes
    for(int i=0;i<nLen;i++) // indexes sample size
    {
        List tempListX1 = basisList[i]; // step into design matrices X of a particular sample size
        int mLen = tempListX1.length(); // number of models
        List returnListSub1(mLen); // part of return list that indexes the generating model
        List tempListY1 = yList[i]; // step into observations of a particular sample size
        for(int j=0;j<mLen;j++) // indexes true generating model
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
                for(int m=0;m<mLen;m++) // indexes regression model to fit for BIC
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
                    for(int d=0;d<noSimDraws;d++) // indexes data generated from a specific regression fit
                    {
                        arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
                        arma::vec tempMean = tempHatMat*tempYs; // fit the current model to the data
                        double tempVar = (accu(pow(tempYs-tempMean,2)))/noObsPts; // estimate the population sd under MLE
                        fillVec[d] = normBICComp(tempYs,tempMean,tempVar,p); // compute the BIC for this fit
                    }
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


