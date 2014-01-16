#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <math.h>

// for shuffling
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

//for distributions
#include "gsl/gsl_cdf.h"

using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b) )

const int experimentCount = 15;
const double smoothingEps = 0.00000000001;
const double trainPart = 0.5;

vector<bool> loadResponse(string fileName)
{
    vector<bool> ret;
    ifstream file;
    file.open(fileName.c_str(),ifstream::in);
    if (!file.good())
    {
        cout<<"error"<<endl;
    }
    cout<<"Loading response"<<endl;

    string word;
    bool temp;
    while(!file.eof())
    {
        file >> temp;
        ret.push_back(temp);
    }
    ret.pop_back();
    file.close();
    cout<<"Loading response finished"<<endl;

    return ret;
}

void loadSamples(vector< vector <double> >& healphy,
                 vector< vector <double> >& cancer,
                 vector<bool>& resp,
                 string fileName)
{
    vector<double> healphyVec;
    vector<double> cancerVec;
    ifstream file;
    file.open(fileName.c_str(),ifstream::in);
    if (!file.good())
    {
        cout<<"error"<<endl;
    }
    cout<<"Loading samples"<<endl;
    
    double temp;
    string line;
    getline(file, line);
    istringstream iss(line);
    int count = 0;
    
    while(!iss.eof())
    {
        iss >> temp;
        if(resp[count])
            cancerVec.push_back(temp);
        else
            healphyVec.push_back(temp);
        count++;
    }
    
    cancer.push_back(cancerVec);
    healphy.push_back(healphyVec);

    while (getline(file, line))
    {
        healphyVec.clear();
        cancerVec.clear();
        istringstream iss(line);
        count = 0;
        while(!iss.eof())
        {
            iss >> temp;
            if(resp[count])
                cancerVec.push_back(temp);
            else
                healphyVec.push_back(temp);
            count++;
        }
        cancer.push_back(cancerVec);
        healphy.push_back(healphyVec);
        
        if(cancer.size() % 10000 == 0)
        {
            cout <<"progress: " << cancer.size()<<endl;
        }
    }
    
    file.close();
    cout<<"Loading samples finished"<<endl;
}

vector<double> mean(vector<vector<double> >& sample, vector<int>& ind)
{
    vector<double> ret;
    double sum = 0;
    vector<double> temp;
    for (int i = 0; i < sample.size(); i++)
    {
        sum = 0;
        for (int j = 0; j < ind.size(); j++)
        {
            sum += sample[i][ind[j]];
        }
        ret.push_back(sum/ind.size());
    }
    return ret;
}

vector<double> sd(vector<vector<double> >& sample, vector<int>& ind)
{
    vector<double> ret;
    double sum = 0;
    double sum2 = 0;
    int N;
    for (int i = 0; i < sample.size(); i++)
    {
        sum = 0;
        sum2 = 0;
        N = ind.size();
        for (int j = 0; j < N; j++)
        {
            sum += sample[i][ind[j]];
            sum2 += sample[i][ind[j]]*sample[i][ind[j]];
        }

        // ret.push_back(sqrt((sum2 - sum*sum/N)/(N-1) + 0.00000000001));
        ret.push_back(sqrt((sum2 - sum*sum/N)/(N) + 0.00000000001));
    }
    return ret;
}

vector<int> get_sample(gsl_rng* randGenerator, int max)
{
    int* vec = new int[max];
    for (int i = 0; i < max; i++)
    {
        vec[i] = i;
    }
    gsl_ran_shuffle(randGenerator, (void*) vec, max,sizeof(*vec));
    vector<int> vecRet;
    for (int i = 0; i < max; i++)
    {
        vecRet.push_back(vec[i]);
    }
    delete [] vec;
    return vecRet;
}

struct SumErrorsStruct
{
    int healphyErr;
    int cancerErr;
    int healphyErrSquared;
    int cancerErrSquared;
    int totalErrSquared;
    int healphyErrTrain;
    int cancerErrTrain;
    SumErrorsStruct():
        healphyErr(0), cancerErr(0), healphyErrSquared(0),
        cancerErrSquared(0), totalErrSquared(0),
        healphyErrTrain(0), cancerErrTrain(0)
        {}
};

struct TrainModel
{
    vector<double> healphyMean;
    vector<double> healphySd;
    vector<double> cancerMean;
    vector<double> cancerSd;
    TrainModel(vector<double> _healphyMean,
               vector<double> _healphySd,
               vector<double> _cancerMean,
               vector<double> _cancerSd):
        healphyMean(_healphyMean),
        healphySd(_healphySd),
        cancerMean(_cancerMean),
        cancerSd(_cancerSd)
        {}
};

struct ResultStatistics
{
    double healphy_mean;
    double healphy_mean_train;
    double healphy_mean_sd;
    double cancer_mean;
    double cancer_mean_train;
    double cancer_mean_sd;
    int tc;
    double total_mean;
    double total_mean_sd;
    string cancerName;
    int experimentCount;
    int healphyTestSize;
    int cancerTestSize;
    void Export(string filename)
        {
            ofstream ofile;
            ofile.open(filename.c_str(), iostream::app);
            ofile<<cancerName<<" ";
            ofile<<experimentCount<<" ";
            ofile<<healphyTestSize<<" ";
            ofile<<healphy_mean<<" ";
            ofile<<healphy_mean_sd<<" ";
            ofile<<healphy_mean_train<<" ";
            ofile<<cancerTestSize<<" ";
            ofile<<cancer_mean<<" ";
            ofile<<cancer_mean_sd<<" ";
            ofile<<cancer_mean_train<<" ";
            ofile<<tc<<" ";
            ofile<<total_mean<<" ";
            ofile<<total_mean_sd<<" "<<endl;
            ofile.close();
        }
};
    

double GetProb(vector<vector<double> >& samples, vector<double>& means, vector<double>& sds,int idx, double subclassProb)
{
    double prob = 0;
    for (int j = 0; j < samples.size(); j++)
    {
        prob += log( gsl_ran_gaussian_pdf(samples[j][idx]-means[j],sds[j])  +smoothingEps);
    }
    prob += subclassProb;
    return prob;
}

void ProcessCancer(string cancerName,gsl_rng* randGenerator, string outFile)
{
    cout<<"Processing "<< cancerName<<" cancer"<<endl;
    // Loading data
    string filenameResp = cancerName + ".resp.csv";
    vector<bool> resp = loadResponse(filenameResp);
    cout<<"resp size = "<<resp.size()<<endl;

    string filenameSample = cancerName + ".sample.csv";
    vector<vector<double> > healphySample,
        cancerSample;
    loadSamples(healphySample, cancerSample, resp, filenameSample);
    cout<<"healphy sample size = "<<healphySample.size()<<endl;
    cout<<"cancer sample size = "<<cancerSample.size()<<endl;

    SumErrorsStruct sumErrs;
    int healphySize = healphySample[0].size(),
        cancerSize = cancerSample[0].size();
    int healphyTrainSize = healphySize*trainPart,
        cancerTrainSize = cancerSize*trainPart,
        healphyTestSize = healphySize - healphyTrainSize,
        cancerTestSize = cancerSize - cancerTrainSize;
    
    // Main loop
    for (int k = 0; k < experimentCount; k++)
    {
        cout<<"##################"<<" experiment index: "<<k<<endl;
        // Form train and test sets
        vector<int > healphyShuffle = get_sample(randGenerator, healphySize),
            cancerShuffle = get_sample(randGenerator, cancerSize);
        vector<int> healphyTrain(healphyShuffle.begin(),
                                 healphyShuffle.begin()+healphyTrainSize),
            healphyTest(healphyShuffle.begin()+healphyTrainSize,
                        healphyShuffle.end()),
            cancerTrain(cancerShuffle.begin(),
                        cancerShuffle.begin()+cancerTrainSize),
            cancerTest(cancerShuffle.begin()+cancerTrainSize,
                       cancerShuffle.end());
        
        TrainModel model(
            mean(healphySample,healphyTrain),
            sd(healphySample,healphyTrain),
            mean(cancerSample,cancerTrain),
            sd(cancerSample,cancerTrain)
            );
        
        // Train complete at this point
        // Predicting
        time_t t = clock();
        int healphyErrLocal = 0;
        int cancerErrLocal = 0;
        double prob0 = 0,
            prob1 = 0,
            healphyClassProb = log(double(healphyTrainSize)/(healphyTrainSize+cancerTrainSize)+smoothingEps),
            cancerClassProb = log(double(cancerTrainSize)/(healphyTrainSize+cancerTrainSize)+smoothingEps);
        cout<<"healphy predicting"<<endl;
        for (int i = 0; i < healphyTestSize; i++)
        {
            prob0 = GetProb(healphySample,model.healphyMean, model.healphySd,healphyTest[i],healphyClassProb);
            prob1 = GetProb(healphySample,model.cancerMean, model.cancerSd,healphyTest[i],cancerClassProb);
//            cout<<prob1- prob0<<endl;
            // cout<<prob1<<endl;
            // cout<<prob0<<endl;

            healphyErrLocal += prob0 < prob1;
        }
        sumErrs.healphyErr += healphyErrLocal;
        sumErrs.healphyErrSquared += healphyErrLocal * healphyErrLocal;

        cout<<"healphy train predicting"<<endl;
        healphyErrLocal = 0;
        for (int i = 0; i < healphyTrainSize; i++)
        {
            prob0 = GetProb(healphySample,model.healphyMean, model.healphySd,healphyTrain[i],healphyClassProb);
            prob1 = GetProb(healphySample,model.cancerMean, model.cancerSd,healphyTrain[i],cancerClassProb);
//            cout<<prob1- prob0<<endl;
            // cout<<prob1<<endl;
            // cout<<prob0<<endl;

            healphyErrLocal += prob0 < prob1;
        }
        sumErrs.healphyErrTrain += healphyErrLocal;
        cout<<"healphy predicting time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
        cout<<"sumErrs = "<<healphyErrLocal<<" out of "<<healphyTestSize<<endl;

        t = clock();
        cancerErrLocal = 0;
        cout<<"cancer predicting"<<endl;
        for (int i = 0; i < cancerTestSize; i++)
        {
            prob0 = GetProb(cancerSample,model.healphyMean, model.healphySd,cancerTest[i],healphyClassProb);
            prob1 = GetProb(cancerSample,model.cancerMean, model.cancerSd,cancerTest[i],cancerClassProb);
            // cout<<prob0 - prob1<<endl;
            cancerErrLocal += prob1 < prob0;
        }
        sumErrs.cancerErr += cancerErrLocal;
        sumErrs.cancerErrSquared += cancerErrLocal * cancerErrLocal;
        
        cout<<"cancer train predicting"<<endl;
        cancerErrLocal = 0;
        for (int i = 0; i < cancerTrainSize; i++)
        {
            prob0 = GetProb(cancerSample,model.healphyMean, model.healphySd,cancerTrain[i],healphyClassProb);
            prob1 = GetProb(cancerSample,model.cancerMean, model.cancerSd,cancerTrain[i],cancerClassProb);
//            cout<<prob1- prob0<<endl;
            // cout<<prob1<<endl;
            // cout<<prob0<<endl;

            cancerErrLocal += prob0 < prob1;
        }
        sumErrs.cancerErrTrain += cancerErrLocal;

        cout<<"cancer predicting time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
        cout<<"sumErrs = "<<cancerErrLocal<<" out of "<<cancerTestSize<<endl;
        sumErrs.totalErrSquared += (cancerErrLocal + healphyErrLocal)*(cancerErrLocal + healphyErrLocal);
    }
    
    cout<<"healphy errors:  "<<sumErrs.healphyErr<< " out of "<<healphyTestSize*experimentCount<<". "<<double(sumErrs.healphyErr)/(healphyTestSize*experimentCount)<<endl;
    cout<<"cancer errors:  "<<sumErrs.cancerErr<< " out of "<<cancerTestSize*experimentCount<<". "<<double(sumErrs.cancerErr)/(cancerTestSize*experimentCount)<<endl;

    // Calculating result statistics and saving it into ResultStatistics structure
    ResultStatistics resStat;
    
    resStat.healphy_mean = double(sumErrs.healphyErr)/(healphyTestSize*experimentCount);
    resStat.healphy_mean_train = double(sumErrs.healphyErrTrain)/(healphyTrainSize*experimentCount);
    resStat.healphy_mean_sd = sqrt(
        double(sumErrs.healphyErrSquared)/(healphyTestSize*
                                           healphyTestSize*
                                           experimentCount) -
        resStat.healphy_mean*resStat.healphy_mean);
    resStat.cancer_mean = double(sumErrs.cancerErr)/(cancerTestSize*experimentCount);
    resStat.cancer_mean_train = double(sumErrs.cancerErrTrain)/(cancerTrainSize*experimentCount);
    resStat.cancer_mean_sd = sqrt(
        double(sumErrs.cancerErrSquared)/(cancerTestSize*
                                          cancerTestSize*
                                          experimentCount) -
        resStat.cancer_mean*resStat.cancer_mean);
    resStat.tc= cancerTestSize + healphyTestSize;
    resStat.total_mean = double(sumErrs.cancerErr+sumErrs.healphyErr)/(resStat.tc*experimentCount);
    resStat.total_mean_sd = sqrt(
        double(sumErrs.totalErrSquared)/(resStat.tc*resStat.tc*
                                         experimentCount) -
        resStat.total_mean*resStat.total_mean);
    resStat.cancerName = cancerName;
    resStat.experimentCount = experimentCount;
    resStat.healphyTestSize = healphyTestSize;
    resStat.cancerTestSize = cancerTestSize;
    resStat.Export(outFile);
}

int main()
{
    gsl_rng* randGenerator= gsl_rng_alloc(gsl_rng_taus);
//     gsl_rng_set(randGenerator, time(NULL));
    gsl_rng_set(randGenerator, 0);
    
    string outFileSuffix("test_train_errs.txt");

    vector<string> cancers;
    cancers.push_back("BLCA");
    // cancers.push_back("READ");
    // cancers.push_back("KIRP");
    // cancers.push_back("LIHC");
    // cancers.push_back("PRAD");
    // cancers.push_back("LUSC");
    // cancers.push_back("COAD");
    // cancers.push_back("LUAD");
    // cancers.push_back("HNSC");
    // cancers.push_back("THCA");
    // cancers.push_back("UCEC");
    // cancers.push_back("KIRC");
    // cancers.push_back("BRCA");
    for (int i = 0; i < cancers.size(); i++)
    {
//        ProcessCancer(cancers[i], randGenerator,cancers[i]+outFileSuffix);
        ProcessCancer(cancers[i], randGenerator, outFileSuffix);
    }
    cout<<"after experiments\n";
    
    gsl_rng_free(randGenerator);
    return 0;
}

// void ProcessCancer1(string cancerName,gsl_rng* randGenerator, string outFile)
// {
//     string filenameResp = cancerName + ".resp.csv";
//     vector<bool> resp = loadResponse(filenameResp);
//     cout<<"resp size = "<<resp.size()<<endl;

//     string filenameSample = cancerName + ".sample.csv";
//     vector<vector<double> > healphySample,
//         cancerSample;
//     loadSamples(healphySample, cancerSample, resp, filenameSample);
//     cout<<"healphy sample size = "<<healphySample.size()<<endl;
//     cout<<"cancer sample size = "<<cancerSample.size()<<endl;
    
//     int experimentCount = 15;
//     double eps = 0.00000000001;
//     int err_healphy = 0;
//     int err_cancer = 0;
//     int kvant_count = MIN(1/0.05,healphySample[0].size()/2);
//     int kvant = (1.0001)/kvant_count * healphySample[0].size()/2;
//     vector<int> err_healphy_glob(kvant_count);
//     for (int i = 0; i < err_healphy_glob.size(); i++)
//         err_healphy_glob[i] = 0;
//     vector<int> err_healphy_glob_train(kvant_count);
//     for (int i = 0; i < err_healphy_glob_train.size(); i++)
//         err_healphy_glob_train[i] = 0;

//     cout<<"healphy train size = "<<healphySample[0].size()/2<<endl;

//     for (int k = 0; k < experimentCount; k++)
//     {
//         vector<int > healphyShuffle(get_sample(randGenerator, healphySample[0].size())),
//             cancerShuffle(get_sample(randGenerator, cancerSample[0].size()));
        
//         vector<int> healphyTest(healphyShuffle.begin()+healphySample[0].size()/2 + 1,
//                               healphyShuffle.end());
//         vector<int> cancerTrain(cancerShuffle.begin(),
//                                    cancerShuffle.begin()+cancerSample[0].size()/2);
//         vector<int> cancerTest(cancerShuffle.begin()+cancerSample[0].size()/2 + 1,
//                                   cancerShuffle.end());
        
//         for (int l = 1; l <= kvant_count; l++)
//         {
//             vector<int> healphyTrain(healphyShuffle.begin(),
//                                         healphyShuffle.begin()+l*kvant);

//             cout<<"experiment "<<k<<"healphyTrain size = "<<healphyTrain.size()<<endl;
        
//             vector<double> means_healphy = mean(healphySample,healphyTrain);
//             vector<double> means_cancer = mean(cancerSample,cancerTrain);
    
//             vector<double> sds_healphy = sd(healphySample,healphyTrain);
//             vector<double> sds_cancer = sd(cancerSample,cancerTrain);

// //train complete at this step
// //testing

//             double prob0 = 0,
//                 prob1 = 0;

//             err_healphy = 0;
//             for (int i = 0; i < healphyTest.size(); i++)
//             {
//                 prob0 = 0;
//                 prob1 = 0;
//                 for (int j = 0; j < healphySample.size(); j++)
//                 {
//                     prob0 += log( gsl_ran_gaussian_pdf(healphySample[j][healphyTest[i]]-means_healphy[j],sds_healphy[j])  +eps);
//                     prob1 += log( gsl_ran_gaussian_pdf(healphySample[j][healphyTest[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
//                 }
//                 err_healphy += prob0 < prob1;
//                 cout<<prob0<<" "<<  prob1<<endl;
//             }
//             cout<<"errs = "<<err_healphy<<" out of "<<healphyTest.size()<<endl;
//             err_healphy_glob[l-1] += err_healphy;
            
//             err_healphy = 0;
//             for (int i = 0; i < healphyTrain.size(); i++)
//             {
//                 prob0 = 0;
//                 prob1 = 0;
//                 for (int j = 0; j < healphySample.size(); j++)
//                 {
//                     prob0 += log( gsl_ran_gaussian_pdf(healphySample[j][healphyTrain[i]]-means_healphy[j],sds_healphy[j])  +eps);
//                     prob1 += log( gsl_ran_gaussian_pdf(healphySample[j][healphyTrain[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
//                 }
//                 err_healphy += prob0 < prob1;
//             }
//             cout<<"errs = "<<err_healphy<<" out of "<<healphyTrain.size()<<endl;
//             err_healphy_glob_train[l-1] += err_healphy;

//         }
//     }
//     vector<double> healphy_mean(kvant_count);
//     for (int i = 0; i < healphy_mean.size(); i++)
//     {
//         cout<<i<<endl;
//         healphy_mean[i] = double(err_healphy_glob[i])/((healphySample[0].size() - healphySample[0].size()/2)*experimentCount);
//     }
    
//     vector<double> healphy_mean_train(kvant_count);
//     for (int i = 0; i < healphy_mean_train.size(); i++)
//     {
//         cout<<i<<endl;
//         healphy_mean_train[i] = double(err_healphy_glob_train[i])/((i+1)*kvant*experimentCount);
//     }
    
//     ofstream ofile;
//     ofile.open(outFile.c_str(), iostream::out);
//     ofile<<kvant<<endl;
//     for (int i = 0; i < healphy_mean.size(); i++)
//     {
//         cout<<i<<endl;
//         ofile<<healphy_mean[i]<<" ";
//     }
//     ofile<<endl;
//     for (int i = 0; i < healphy_mean_train.size(); i++)
//     {
//         cout<<i<<endl;
//         ofile<<healphy_mean_train[i]<<" ";
//     }
//     ofile<<endl;
    
//     ofile.close();
//     cout<<"before closing\n";
// }
