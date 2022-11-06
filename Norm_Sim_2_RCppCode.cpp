#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
#define pi 3.14159265358979
// [[Rcpp::depends(RcppArmadillo)]]


// code to simulate normal data
// [[Rcpp::export]]
List normDatSim(arma::vec ns,List MLEList,int noSimDraws)
{
    int nLen = MLEList.length(); // get the number of different n's
    List returnList(nLen); // the return list that indexes the sample sizes
    for(int i=0;i<nLen;++i) // index the sample sizes
    {
        int tempN = ns(i);
        List MLEsOfFixedN = MLEList[i];
        int mLen=MLEsOfFixedN.length(); // get the number of different models
        List returnListSub1(mLen); // part of return list that indexes the generating model
        for(int j=0;j<mLen;++j)
        {
            arma::mat MLEsOfFixedNAndModel = MLEsOfFixedN[j];
            int noDraws =  MLEsOfFixedNAndModel.n_cols;
            List returnListSub2(noDraws); // part of return list that indexes the exact fits for each model
            for(int k=0;k<noDraws;++k)
            {
                arma::vec tempMLEs = MLEsOfFixedNAndModel.col(k);
                double tempMean = tempMLEs[0];
                double tempSig = tempMLEs[1];
                arma::mat tempNormDraws(tempN,noSimDraws,fill::randn);
                tempNormDraws=tempNormDraws*tempSig+tempMean; // scale and translate values
                returnListSub2[k] = tempNormDraws;
            }
            returnListSub1[j]=returnListSub2;
        }
        returnList[i]=returnListSub1;
        Rcout << "Iteration " << i+1 << " Complete" << endl;
    }

    return(returnList);
}

// compute the BIC for a normally distributed random variable
// [[Rcpp::export]]
double normBICComp(arma::vec ys, string model)
{
    double N = ys.n_elem;
    double BIC;
    if(model=="N(0,1)")
    {
        arma::vec sqrdDiffs = pow(ys-0,2);
        double LL = -N/2.0*log(2.0*pi)-N/2.0*log(1.0)-accu(sqrdDiffs)/(2.0*1.0);
        BIC = -2.0*LL;
    }
    else if(model=="N(mu,1)")
    {
        double MLEMean = mean(ys);
        arma::vec sqrdDiffs = pow(ys-MLEMean,2);
        double LL = -N/2.0*log(2.0*pi)-N/2.0*log(1.0)-accu(sqrdDiffs)/(2.0*1.0);
        BIC = -2.0*LL + 1*log(N);
    }
    else // for the N(mu,sig) model
    {
        double MLEMean = mean(ys);
        arma::vec sqrdDiffs = pow(ys-MLEMean,2);
        double accuSqrdDiffs = accu(sqrdDiffs);
        double MLEVar = accuSqrdDiffs/N;
        double LL = -N/2.0*log(2.0*pi)-N/2.0*log(MLEVar)-accuSqrdDiffs/(2.0*MLEVar);
        BIC = -2.0*LL + 2.0*log(N);
    }
    return(BIC);
}

// compute the BIC for a normally distributed random variable
// [[Rcpp::export]]
arma::vec normBICComp2(arma::mat tempYMat,string model,arma::vec fillVec, int noSimDraws)
{
    double N = tempYMat.n_rows;
    if(model=="N(0,1)")
    {
        for(int d=0;d<noSimDraws;d++) // indexes data generated from a specific regression fit
        {
            arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
            arma::vec sqrdDiffs = pow(tempYs-0,2);
            double LL = -N/2.0*log(2.0*pi)-N/2.0*log(1.0)-accu(sqrdDiffs)/(2.0*1.0);
            fillVec[d] = -2.0*LL; // compute the BIC for this fit
        }
    }
    else if(model=="N(mu,1)")
    {
        for(int d=0;d<noSimDraws;d++) // indexes data generated from a specific regression fit
        {
            arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
            double MLEMean = mean(tempYs);
            arma::vec sqrdDiffs = pow(tempYs-MLEMean,2);
            double LL = -N/2.0*log(2.0*pi)-N/2.0*log(1.0)-accu(sqrdDiffs)/(2.0*1.0);
            fillVec[d] = -2.0*LL + 1.0*log(N); // compute the BIC for this fit
        }
    }
    else // for the N(mu,sig) model
    {
        for(int d=0;d<noSimDraws;d++) // indexes data generated from a specific regression fit
        {
            arma::vec tempYs = tempYMat.col(d); // pull a particular simulated sample
            double MLEMean = mean(tempYs);
            arma::vec sqrdDiffs = pow(tempYs-MLEMean,2);
            double accuSqrdDiffs = accu(sqrdDiffs);
            double MLEVar = accuSqrdDiffs/N;
            double LL = -N/2.0*log(2.0*pi)-N/2.0*log(MLEVar)-accuSqrdDiffs/(2.0*MLEVar);
            fillVec[d] =  -2.0*LL + 2.0*log(N); // compute the BIC for this fit
        }
    }
    return(fillVec);
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


// [[Rcpp::export]]
List RcppICComps(List yList,Rcpp::StringVector models)
{
    int nLen = yList.length(); // number of sample sizes
    List returnList(nLen); // the return list that indexes the sample sizes
    for(int i=0;i<nLen;i++) // indexes sample size
    {
        int mLen = models.length(); // number of models
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
                for(int m=0;m<mLen;m++) // indexes normal model to fit for BIC
                {
                    String curModel = models[m];
                    // vector to fill the BIC computation for each simulated draw generated from model j but assuming here that model m is true
                    arma::vec fillVec = vec(noSimDraws);
                    fillVec=normBICComp2(tempYMat,curModel,fillVec,noSimDraws);
                    fillMat.col(m) = fillVec; // column m of fillMat will hold the BICs for each draw computed under model m
                }
                returnListSub2[k] = fillMat; // store all BIC computations for normal fit k
            }
        returnListSub1[j] = returnListSub2; // store the list of BIC computations for all noDraws matrices truly generated from model j
        }
    returnList[i] = returnListSub1; // store the list of BIC computations for all noDraws matrices truly generated from model j with sample sizes i
    Rcout << "Its. for Sample Size Index " << i+1 << " completed" << endl;
    }
    return(returnList);
}

// [[Rcpp::export]]
List RcppSimDGOFs(List yList,Rcpp::StringVector models)
{
    int nLen = yList.length(); // number of sample sizes
    List returnList(nLen); // the return list that indexes the sample sizes
    for(int i=0;i<nLen;i++) // indexes sample size
    {
        int mLen = models.length(); // number of models
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
                for(int m=0;m<mLen;m++) // indexes normal model to fit for BIC
                {
                    String curModel = models[m];
                    // vector to fill the BIC computation for each simulated draw generated from model j but assuming here that model m is true
                    arma::vec fillVec = vec(noSimDraws);
                    fillVec=normBICComp2(tempYMat,curModel,fillVec,noSimDraws);
                    fillMat.col(m) = fillVec; // column m of fillMat will hold the BICs for each draw computed under model m
                }
                arma::mat tempDGOFs = DGOFCompforMatrix(fillMat);
                returnListSub2[k] = tempDGOFs; // store all BIC computations for normal fit k
            }
        returnListSub1[j] = returnListSub2; // store the list of BIC computations for all noDraws matrices truly generated from model j
        }
    returnList[i] = returnListSub1; // store the list of BIC computations for all noDraws matrices truly generated from model j with sample sizes i
    Rcout << "Its. for Sample Size Index " << i+1 << " completed" << endl;
    }
    return(returnList);
}