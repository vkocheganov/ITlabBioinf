#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;
#define MIN(a,b) ((a) < (b) ? (a) : (b) )


vector<bool> loadResponse(string fileName)
{
    vector<bool> ret;
    ifstream file;
    file.open(fileName.c_str(),ifstream::in);
    if (!file.good())
    {
        cout<<"error"<<endl;
    }
    
    string word;
    bool temp;
    while(!file.eof())
    {
        file >> temp;
        ret.push_back(temp);
//        cout <<temp<<endl;
    }
    ret.pop_back();
    file.close();
    return ret;
}

void loadSamples(vector< vector <double> >& healphy,
                 vector< vector <double> >& cancer,
                 vector<bool>& resp,
                 string fileName)
{
    healphy.clear();
    cancer.clear();
    vector<double> retHealphy;
    vector<double> retCancer;
    ifstream file;
    int size;
    file.open(fileName.c_str(),ifstream::in);
    if (!file.good())
    {
        cout<<"error"<<endl;
    }
    
    double temp;
    string line;
    retHealphy.clear();
    retCancer.clear();
    getline(file, line);
    istringstream iss(line);
    int count = 0;
    
    while(!iss.eof())
    {
        iss >> temp;
        if(resp[count])
            retCancer.push_back(temp);
        else
            retHealphy.push_back(temp);
        count++;
    }
    
    cancer.push_back(retCancer);
    healphy.push_back(retHealphy);

    // process pair (a,b)

    while (getline(file, line))
    {
        retHealphy.clear();
        retCancer.clear();
        istringstream iss(line);
        count = 0;
        while(!iss.eof())
        {
            iss >> temp;
            if(resp[count])
                retCancer.push_back(temp);
            else
                retHealphy.push_back(temp);
            count++;
        }
        cancer.push_back(retCancer);
        healphy.push_back(retHealphy);
        
        if(cancer.size() % 10000 == 0)
        {
            cout <<cancer.size()<<endl;
        }
        
        // process pair (a,b)
    }
    
    file.close();
}

#include <algorithm>
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
    // for (int i = 0; i < sample.size(); i++)
    // {
    //     sum = 0;
    //     temp.clear();
    //     for (int j = 0; j < ind.size(); j++)
    //     {
    //         temp.push_back(sample[i][ind[j]]);
    //     }
    //     sort(temp.begin(),temp.end());
    //     ret.push_back(temp[temp.size()/2]);
    // }

    return ret;
}

vector<double> var(vector<vector<double> >& sample, vector<int>& ind)
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
        ret.push_back((sum2 - sum*sum/N)/(N));
    }
    return ret;
}

#include <math.h>
vector<double> sd(vector<vector<double> >& sample, vector<int>& ind)
//vector<double> sd(vector<vector<double> >& sample, vector<int>& ind, vector<double>& means)
{
    vector<double> ret;
    double sum = 0;
//    double temp;
    double sum2 = 0;
    int N;
    // for (int i = 0; i < sample.size(); i++)
    // {
    //     sum = 0;
    //     N = ind.size();
    //     for (int j = 0; j < N; j++)
    //     {
    //         temp = means[i] - sample[i][ind[j]];
    //         sum += temp*temp;
    //     }
    //     ret.push_back(sqrt(sum/N));
    // }
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

        ret.push_back(sqrt((sum2 - sum*sum/N)/(N) + 0.00000000001));
    }
    return ret;
}

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h" // for shuffling

vector<int> get_sample(gsl_rng* r, int max)
{
    int* vec = new int[max];
    for (int i = 0; i < max; i++)
    {
        vec[i] = i;
    }
    gsl_ran_shuffle(r, (void*) vec, max,sizeof(*vec));
    vector<int> vecRet;
    for (int i = 0; i < max; i++)
    {
        vecRet.push_back(vec[i]);
    }
    delete [] vec;
    return vecRet;
}

//g++ -O3 bayes.cpp  -lgsl -lgslcblas  -o bayes && ./bayes
#include "gsl/gsl_cdf.h"

void make_experiment(string cancerName,gsl_rng* r, string outFile)
{
    string filenameResp = cancerName + ".resp.csv";
    string filenameSample = cancerName + ".sample.csv";
    
    vector<bool> resp = loadResponse(filenameResp);
    cout<<"resp size = "<<resp.size()<<endl;

    vector<vector<double> > healphy_sample,
        cancer_sample;
    loadSamples(healphy_sample, cancer_sample, resp, filenameSample);
    cout<<"healphy sample size = "<<healphy_sample.size()<<endl;
    cout<<"cancer sample size = "<<cancer_sample.size()<<endl;
    


    int experimentCount = 15;
    double eps = 0.00000000001;
    int err_healphy = 0;
    int err_cancer = 0;
    int err_healphy_glob = 0;
    int err_cancer_glob = 0;
    int err_healphy_glob_train = 0;
    int err_cancer_glob_train = 0;
    int err_healphy_glob2 = 0;
    int err_cancer_glob2 = 0;
    int err_total_glob2 = 0;
    vector<int > healphyShuffle,
        cancerShuffle;
    vector<int> healphy_train;
    vector<int> healphy_test;
    vector<int> cancer_train;
    vector<int> cancer_test;
    for (int k = 0; k < experimentCount; k++)
    {
        healphyShuffle = get_sample(r, healphy_sample[0].size());
        cancerShuffle = get_sample(r, cancer_sample[0].size());
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
        // vector<double> sds_healphy = sd(healphy_sample,healphy_train,means_healphy);
        // vector<double> sds_cancer = sd(cancer_sample,cancer_train, means_cancer);
        vector<double> sds_healphy = sd(healphy_sample,healphy_train);
        vector<double> sds_cancer = sd(cancer_sample,cancer_train);
        cout<<"sd time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;

//train complete at this step
//testing

        double prob0 = 0,
            prob1 = 0;

        // err_healphy = 0;
        // for (int i = 0; i < healphy_train.size(); i++)
        // {
        //     prob0 = 0;
        //     prob1 = 0;
        //     for (int j = 0; j < healphy_sample.size(); j++)
        //     {
        //         prob0 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_train[i]]-means_healphy[j],sds_healphy[j])  +eps);
        //         prob1 += log( gsl_ran_gaussian_pdf(healphy_sample[j][healphy_train[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
        //     }
        //     prob0 += log(double(healphy_train.size())/(healphy_train.size()+cancer_train.size())+eps);
        //     prob1 += log(double(cancer_train.size())/(healphy_train.size()+cancer_train.size())+eps);
        //     err_healphy += prob0 < prob1;
        // }
        // cout<<"healphy train errs = "<<err_healphy<<" out of "<<healphy_train.size()<<endl;
        // err_healphy_glob_train += err_healphy;
        
        t = clock();
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
            // prob0 += log(double(healphy_train.size())/(healphy_train.size()+cancer_train.size())+eps);
            // prob1 += log(double(cancer_train.size())/(healphy_train.size()+cancer_train.size())+eps);
//            cout<<prob1 - prob0<<endl;
            err_healphy += prob0 < prob1;
        }
        cout<<"healphy predicting time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
        cout<<"errs = "<<err_healphy<<" out of "<<healphy_test.size()<<endl;


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
//         err_cancer_glob_train += err_cancer;

        cout<<"errs cancer train = "<<err_cancer<<" out of "<<cancer_train.size()<<endl;
        t = clock();
        err_cancer = 0;
        for (int i = 0; i < cancer_test.size(); i++)
        {
            prob0 = 0;
            prob1 = 0;
            for (int j = 0; j < cancer_sample.size(); j++)
            {
                prob0 += log( gsl_ran_gaussian_pdf(cancer_sample[j][cancer_test[i]]-means_healphy[j],sds_healphy[j])  +eps);
                prob1 += log( gsl_ran_gaussian_pdf(cancer_sample[j][cancer_test[i]]-means_cancer[j] ,sds_cancer[j]) +eps);
            }
            // prob0 += log(double(healphy_train.size())/(cancer_train.size()+healphy_train.size())+eps);
            // prob1 += log(double(cancer_train.size())/(cancer_train.size()+healphy_train.size())+eps);
//            cout<<prob0 - prob1<<endl;
            err_cancer += prob1 < prob0;
        }
        err_healphy_glob += err_healphy;
        err_healphy_glob2 += err_healphy * err_healphy;
        err_cancer_glob += err_cancer;
        err_cancer_glob2 += err_cancer * err_cancer;
        err_total_glob2 += (err_cancer + err_healphy)*(err_cancer + err_healphy);
        cout<<"cancer predicting time = "<<double(clock() - t)/CLOCKS_PER_SEC<<endl;
        cout<<"errs = "<<err_cancer<<" out of "<<cancer_test.size()<<endl;
    }
    
    cout<<"healphy errors:  "<<err_healphy_glob<< " out of "<<healphy_test.size()*experimentCount<<". "<<double(err_healphy_glob)/(healphy_test.size()*experimentCount)<<endl;
    cout<<"cancer errors:  "<<err_cancer_glob<< " out of "<<cancer_test.size()*experimentCount<<". "<<double(err_cancer_glob)/(cancer_test.size()*experimentCount)<<endl;
    
    double healphy_mean = double(err_healphy_glob)/(healphy_test.size()*experimentCount);
    double healphy_mean_train = double(err_healphy_glob_train)/(healphy_train.size()*experimentCount);
    double healphy_mean_sd = sqrt(
        double(err_healphy_glob2)/(healphy_test.size()*
                                   healphy_test.size()*
                                   experimentCount) -
        healphy_mean*healphy_mean);
    double cancer_mean = double(err_cancer_glob)/(cancer_test.size()*experimentCount);
    double cancer_mean_train = double(err_cancer_glob_train)/(cancer_train.size()*experimentCount);
    double cancer_mean_sd = sqrt(
        double(err_cancer_glob2)/(cancer_test.size()*
                                   cancer_test.size()*
                                   experimentCount) -
        cancer_mean*cancer_mean);
    int tc= cancer_test.size() + healphy_test.size();
    double total_mean = double(err_cancer_glob+err_healphy_glob)/(tc*experimentCount);
    double total_mean_sd = sqrt(
        double(err_total_glob2)/(tc*tc*
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

void make_experiment1(string cancerName,gsl_rng* rand_generator, string outFile)
{
    string filenameResp = cancerName + ".resp.csv";
    string filenameSample = cancerName + ".sample.csv";
    
    vector<bool> resp = loadResponse(filenameResp);
    cout<<"resp size = "<<resp.size()<<endl;

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
    gsl_rng_set(rand_generator, time(NULL));
    
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
        make_experiment1(cancers[i], rand_generator,cancers[i]+outFileSuffix);
    }
    cout<<"after experiments\n";
    
    gsl_rng_free(rand_generator);
    return 0;
}


