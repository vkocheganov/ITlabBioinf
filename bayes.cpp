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

vector<int> get_sample(gsl_rng* rand_generator, int max)
{
    int* vec = new int[max];
    for (int i = 0; i < max; i++)
    {
        vec[i] = i;
    }
    gsl_ran_shuffle(rand_generator, (void*) vec, max,sizeof(*vec));
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
    
void ProcessCancer(string cancerName,gsl_rng* rand_generator, string outFile)
{
    // Loading data
    string filenameResp = cancerName + ".resp.csv";
    vector<bool> resp = loadResponse(filenameResp);
    cout<<"resp size = "<<resp.size()<<endl;

    string filenameSample = cancerName + ".sample.csv";
    vector<vector<double> > healphy_sample,
        cancer_sample;
    loadSamples(healphy_sample, cancer_sample, resp, filenameSample);
    cout<<"healphy sample size = "<<healphy_sample.size()<<endl;
    cout<<"cancer sample size = "<<cancer_sample.size()<<endl;

    // Main loop
    SumErrorsStruct errs;
    
    vector<int > healphyShuffle,
        cancerShuffle;
    vector<int> healphy_train;
    vector<int> healphy_test;
    vector<int> cancer_train;
    vector<int> cancer_test;
    for (int k = 0; k < experimentCount; k++)
    {
        int healphyErrLocal = 0;
        int cancerErrLocal = 0;
        healphyShuffle = get_sample(rand_generator, healphy_sample[0].size());
        cancerShuffle = get_sample(rand_generator, cancer_sample[0].size());
        healphy_train = vector<int>(healphyShuffle.begin(),
                               healphyShuffle.begin()+healphy_sample[0].size()/2);
        healphy_test = vector<int>(healphyShuffle.begin()+healphy_sample[0].size()/2 + 1,
                              healphyShuffle.end());

        cancer_train = vector<int>(cancerShuffle.begin(),
                              cancerShuffle.begin()+cancer_sample[0].size()/2);
        cancer_test = vector<int>(cancerShuffle.begin()+cancer_sample[0].size()/2 + 1,
                             cancerShuffle.end());
        cout<<"experiment "<<k<<endl;
        
        time_t t = clock();
        vector<double> means_healphy = mean(healphy_sample,healphy_train);
        vector<double> means_cancer = mean(cancer_sample,cancer_train);
        cout<<"mean time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
    
        t = clock();
        vector<double> sds_healphy = sd(healphy_sample,healphy_train);
        vector<double> sds_cancer = sd(cancer_sample,cancer_train);
        cout<<"sd time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;

//train complete at this step
//testing

        double prob0 = 0,
            prob1 = 0;
        
        t = clock();
        healphyErrLocal = 0;
        for (int i = 0; i < healphy_test.size(); i++)
        {
            prob0 = 0;
            prob1 = 0;
            for (int j = 0; j < healphy_sample.size(); j++)
            {
                prob0 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_test[i]]-means_healphy[j],sds_healphy[j])  +smoothingEps);
                prob1 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_test[i]]-means_cancer[j] ,sds_cancer[j]) +smoothingEps);
            }
            prob0 += log(double(healphy_train.size())/(healphy_train.size()+cancer_train.size())+smoothingEps);
            prob1 += log(double(cancer_train.size())/(healphy_train.size()+cancer_train.size())+smoothingEps);
//            cout<<prob1 - prob0<<endl;
            healphyErrLocal += prob0 < prob1;
        }
        cout<<"healphy predicting time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
        cout<<"errs = "<<healphyErrLocal<<" out of "<<healphy_test.size()<<endl;


//         err_cancer = 0;
//         for (int i = 0; i < cancer_train.size(); i++)
//         {
//             prob0 = 0;
//             prob1 = 0;
//             for (int j = 0; j < cancer_sample.size(); j++)
//             {
//                 prob0 += log( gsl_ran_gaussian_pdf(cancer_sample[j][cancer_train[i]]-means_healphy[j],sds_healphy[j])  +eps);
//                 prob1 += log( gsl_ran_gaussian_pdf(cancer_sample[j][cancer_train[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
//             }
//             prob0 += log(double(healphy_train.size())/(cancer_train.size()+healphy_train.size())+eps);
//             prob1 += log(double(cancer_train.size())/(cancer_train.size()+healphy_train.size())+eps);
// //            cout<<prob0 - prob1<<endl;
//             err_cancer += prob1 < prob0;
//         }
//         errs.cancerErr_train += err_cancer;

        cout<<"errs cancer train = "<<cancerErrLocal<<" out of "<<cancer_train.size()<<endl;
        t = clock();
        cancerErrLocal = 0;
        for (int i = 0; i < cancer_test.size(); i++)
        {
            prob0 = 0;
            prob1 = 0;
            for (int j = 0; j < cancer_sample.size(); j++)
            {
                prob0 += log( gsl_ran_gaussian_pdf(cancer_sample[j][cancer_test[i]]-means_healphy[j],sds_healphy[j])  +smoothingEps);
                prob1 += log( gsl_ran_gaussian_pdf(cancer_sample[j][cancer_test[i]]-means_cancer[j] ,sds_cancer[j]) +smoothingEps);
            }
            // prob0 += log(double(healphy_train.size())/(cancer_train.size()+healphy_train.size())+eps);
            // prob1 += log(double(cancer_train.size())/(cancer_train.size()+healphy_train.size())+eps);
//            cout<<prob0 - prob1<<endl;
            cancerErrLocal += prob1 < prob0;
        }
        errs.healphyErr += healphyErrLocal;
        errs.healphyErrSquared += healphyErrLocal * healphyErrLocal;
        errs.cancerErr += cancerErrLocal;
        errs.cancerErrSquared += cancerErrLocal * cancerErrLocal;
        errs.totalErrSquared += (cancerErrLocal + healphyErrLocal)*(cancerErrLocal + healphyErrLocal);
        cout<<"cancer predicting time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
        cout<<"errs = "<<cancerErrLocal<<" out of "<<cancer_test.size()<<endl;
    }
    
    cout<<"healphy errors:  "<<errs.healphyErr<< " out of "<<healphy_test.size()*experimentCount<<". "<<double(errs.healphyErr)/(healphy_test.size()*experimentCount)<<endl;
    cout<<"cancer errors:  "<<errs.cancerErr<< " out of "<<cancer_test.size()*experimentCount<<". "<<double(errs.cancerErr)/(cancer_test.size()*experimentCount)<<endl;
    
    double healphy_mean = double(errs.healphyErr)/(healphy_test.size()*experimentCount);
    double healphy_mean_train = double(errs.healphyErrTrain)/(healphy_train.size()*experimentCount);
    double healphy_mean_sd = sqrt(
        double(errs.healphyErrSquared)/(healphy_test.size()*
                                   healphy_test.size()*
                                   experimentCount) -
        healphy_mean*healphy_mean);
    double cancer_mean = double(errs.cancerErr)/(cancer_test.size()*experimentCount);
    double cancer_mean_train = double(errs.cancerErrTrain)/(cancer_train.size()*experimentCount);
    double cancer_mean_sd = sqrt(
        double(errs.cancerErrSquared)/(cancer_test.size()*
                                   cancer_test.size()*
                                   experimentCount) -
        cancer_mean*cancer_mean);
    int tc= cancer_test.size() + healphy_test.size();
    double total_mean = double(errs.cancerErr+errs.healphyErr)/(tc*experimentCount);
    double total_mean_sd = sqrt(
        double(errs.totalErrSquared)/(tc*tc*
                                  experimentCount) -
        total_mean*total_mean);
    
    ofstream ofile;
    ofile.open(outFile.c_str(), iostream::app);
    ofile<<cancerName<<" ";
    ofile<<experimentCount<<" ";
    ofile<<healphy_test.size()<<" ";
    ofile<<healphy_mean<<" ";
    ofile<<healphy_mean_sd<<" ";
    ofile<<healphy_mean_train<<" ";
    ofile<<cancer_test.size()<<" ";
    ofile<<cancer_mean<<" ";
    ofile<<cancer_mean_sd<<" ";
    ofile<<cancer_mean_train<<" ";
    ofile<<tc<<" ";
    ofile<<total_mean<<" ";
    ofile<<total_mean_sd<<" "<<endl;
    ofile.close();
}

void ProcessCancer1(string cancerName,gsl_rng* rand_generator, string outFile)
{
    string filenameResp = cancerName + ".resp.csv";
    vector<bool> resp = loadResponse(filenameResp);
    cout<<"resp size = "<<resp.size()<<endl;

    string filenameSample = cancerName + ".sample.csv";
    vector<vector<double> > healphy_sample,
        cancer_sample;
    loadSamples(healphy_sample, cancer_sample, resp, filenameSample);
    cout<<"healphy sample size = "<<healphy_sample.size()<<endl;
    cout<<"cancer sample size = "<<cancer_sample.size()<<endl;
    
    int experimentCount = 15;
    double eps = 0.00000000001;
    int err_healphy = 0;
    int err_cancer = 0;
    int kvant_count = MIN(1/0.05,healphy_sample[0].size()/2);
    int kvant = (1.0001)/kvant_count * healphy_sample[0].size()/2;
    vector<int> err_healphy_glob(kvant_count);
    for (int i = 0; i < err_healphy_glob.size(); i++)
        err_healphy_glob[i] = 0;
    vector<int> err_healphy_glob_train(kvant_count);
    for (int i = 0; i < err_healphy_glob_train.size(); i++)
        err_healphy_glob_train[i] = 0;

    cout<<"healphy train size = "<<healphy_sample[0].size()/2<<endl;

    for (int k = 0; k < experimentCount; k++)
    {
        vector<int > healphyShuffle(get_sample(rand_generator, healphy_sample[0].size())),
            cancerShuffle(get_sample(rand_generator, cancer_sample[0].size()));
        
        vector<int> healphy_test(healphyShuffle.begin()+healphy_sample[0].size()/2 + 1,
                              healphyShuffle.end());
        vector<int> cancer_train(cancerShuffle.begin(),
                                   cancerShuffle.begin()+cancer_sample[0].size()/2);
        vector<int> cancer_test(cancerShuffle.begin()+cancer_sample[0].size()/2 + 1,
                                  cancerShuffle.end());
        
        for (int l = 1; l <= kvant_count; l++)
        {
            vector<int> healphy_train(healphyShuffle.begin(),
                                        healphyShuffle.begin()+l*kvant);

            cout<<"experiment "<<k<<"healphy_train size = "<<healphy_train.size()<<endl;
        
            vector<double> means_healphy = mean(healphy_sample,healphy_train);
            vector<double> means_cancer = mean(cancer_sample,cancer_train);
    
            vector<double> sds_healphy = sd(healphy_sample,healphy_train);
            vector<double> sds_cancer = sd(cancer_sample,cancer_train);

//train complete at this step
//testing

            double prob0 = 0,
                prob1 = 0;

            err_healphy = 0;
            for (int i = 0; i < healphy_test.size(); i++)
            {
                prob0 = 0;
                prob1 = 0;
                for (int j = 0; j < healphy_sample.size(); j++)
                {
                    prob0 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_test[i]]-means_healphy[j],sds_healphy[j])  +eps);
                    prob1 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_test[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
                }
                err_healphy += prob0 < prob1;
                cout<<prob0<<" "<<  prob1<<endl;
            }
            cout<<"errs = "<<err_healphy<<" out of "<<healphy_test.size()<<endl;
            err_healphy_glob[l-1] += err_healphy;
            
            err_healphy = 0;
            for (int i = 0; i < healphy_train.size(); i++)
            {
                prob0 = 0;
                prob1 = 0;
                for (int j = 0; j < healphy_sample.size(); j++)
                {
                    prob0 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_train[i]]-means_healphy[j],sds_healphy[j])  +eps);
                    prob1 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_train[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
                }
                err_healphy += prob0 < prob1;
            }
            cout<<"errs = "<<err_healphy<<" out of "<<healphy_train.size()<<endl;
            err_healphy_glob_train[l-1] += err_healphy;

        }
    }
    vector<double> healphy_mean(kvant_count);
    for (int i = 0; i < healphy_mean.size(); i++)
    {
        cout<<i<<endl;
        healphy_mean[i] = double(err_healphy_glob[i])/((healphy_sample[0].size() - healphy_sample[0].size()/2)*experimentCount);
    }
    
    vector<double> healphy_mean_train(kvant_count);
    for (int i = 0; i < healphy_mean_train.size(); i++)
    {
        cout<<i<<endl;
        healphy_mean_train[i] = double(err_healphy_glob_train[i])/((i+1)*kvant*experimentCount);
    }
    
    ofstream ofile;
    ofile.open(outFile.c_str(), iostream::out);
    ofile<<kvant<<endl;
    for (int i = 0; i < healphy_mean.size(); i++)
    {
        cout<<i<<endl;
        ofile<<healphy_mean[i]<<" ";
    }
    ofile<<endl;
    for (int i = 0; i < healphy_mean_train.size(); i++)
    {
        cout<<i<<endl;
        ofile<<healphy_mean_train[i]<<" ";
    }
    ofile<<endl;
    
    ofile.close();
    cout<<"before closing\n";
}

int main()
{
    gsl_rng* rand_generator= gsl_rng_alloc(gsl_rng_taus);
    // gsl_rng_set(rand_generator, time(NULL));
    gsl_rng_set(rand_generator, 0);
    
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
        ProcessCancer(cancers[i], rand_generator,cancers[i]+outFileSuffix);
    }
    cout<<"after experiments\n";
    
    gsl_rng_free(rand_generator);
    return 0;
}

//g++ -O3 bayes.cpp  -lgsl -lgslcblas  -o bayes && ./bayes
