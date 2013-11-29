#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

vector<double> loadResponse(string fileName)
{
    vector<double> ret;
    ifstream file;
    file.open(fileName.c_str(),ifstream::in);
    if (!file.good())
    {
        cout<<"error"<<endl;
    }
    
    string word;
    double temp;
    while(!file.eof())
    {
        file >> temp;
        ret.push_back(temp);
        cout <<temp<<endl;
    }
    ret.pop_back();
    file.close();
    return ret;
}

vector< vector <double> > loadSamples(string fileName)
{
    vector<vector<double> > ret;
    vector<double> ret1;
    ifstream file;
    int size;
    file.open(fileName.c_str(),ifstream::in);
    if (!file.good())
    {
        cout<<"error"<<endl;
    }
    
    double temp;
    string line;
    ret1.clear();
    getline(file, line);
    istringstream iss(line);
    while(!iss.eof())
    {
        iss >> temp;
//        cout<< temp<<" ";
        ret1.push_back(temp);
    }
//    cout<<endl;
    ret.push_back(ret1);
    size = ret1.size();
    // process pair (a,b)

    while (getline(file, line))
    {
        ret1.clear();
        istringstream iss(line);
        while(!iss.eof())
        {
            iss >> temp;
            cout<< temp<<" ";
            ret1.push_back(temp);
        }
        cout<<endl;
        ret.push_back(ret1);
        if (size != ret1.size())
        {
            cout<<"!!!!!!!!!!!!!!!!!!!!!! size mismatched !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            return ret;
        }
        if(ret.size() % 10000 == 0)
            cout <<ret.size()<<endl;
        // process pair (a,b)
    }
    
    file.close();
    cout <<"sample size = "<< size<<endl;
    return ret;
}


int main()
{
    cout<<"hello\n";
//    vector<double> v = loadResponse("/home/maths/Documents/betta_data/ITlabBioinf/testResp.csv");
    vector<double> v = loadResponse("../resp.csv");
    cout<<"size = "<<v.size()<<endl;
    loadSamples("../data.csv");
//    loadSamples("/home/maths/Documents/betta_data/ITlabBioinf/testSamples.csv");
//    cout<<"size = "<<v.size()<<endl;
    return 0;
}
