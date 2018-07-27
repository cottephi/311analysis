#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <TMyFileHeader.h>


using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}


void GetDataFromRun(vector<int> run_list = {}){
  
  gErrorIgnoreLevel = kError;

  string path_header = runs_headers;
  string path_slow_control = slow_control;
  if (run_list.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path_header << "*..." << endl;
    #endif
    string wildcard_path = path_header + "*";
    for( auto irun : glob(wildcard_path) ){
      run_list.push_back(atoi(irun.data()));
    }
  }
  
  for(auto run : run_list){
    string file_header = path_header + to_string(run) + ".root";
    string directory = path_slow_control + to_string(run) + "/";
    if(!ExistTest(string(directory).data())){
      #if verbose
      cout << "Creating directory " << directory << "..." << endl;
      #endif
      mkdir(string(directory).data(),0777);
    }
    string file_slow_control = directory + to_string(run) + ".txt";
    if(!ExistTest(file_header)){cout << "ERROR: file " << file_header << " not found" << endl; return;}
    
    
    TFile ifile(file_header.data(), "READ");
    TMyFileHeader *header = (TMyFileHeader*)((TKey*)ifile.GetListOfKeys()->At(0))->ReadObj();
    int tstart = header->GetStartTime();
    int tend = header->GetEndTime();
    string command  = path_311analysis+"slow_control/Pcotte_getdata.sh " + to_string(tstart) + " " + to_string(tend) + " " + to_string(run);
    system(command.data());
    ifile.Close();
  }
}
