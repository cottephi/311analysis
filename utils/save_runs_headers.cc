#include <stdlib.h>
#include <stdio.h>
#include <311Lib.h>
#include <TMyFileHeader.h>

using namespace std;

void save_runs_headers(vector<int> run_list = {}, string version = "Feb"){
  if(!Load_Version(version)){return;}
  if(!load_run_lists()){return;}
  
  string path = path_wa105_311data;
  string path_output = runs_headers;
  if (run_list.size() == 0){
    path = path_wa105_311data;
    #if verbose
    cout << "Processing all runs in " << path_wa105_311data << "*..." << endl;
    #endif
    string wildcard_path = path_wa105_311data + "*";
    for( auto irun : glob(wildcard_path) ){
      run_list.push_back(atoi(irun.data()));
    }
  }
  string filename = "";
  
  for(auto run : run_list){
  
    string ofilename = path_output+to_string(run)+".root";
    double drift = runs_and_fields[run]["Drift"];
    double amplification = runs_and_fields[run]["Extraction"];
    double extraction = runs_and_fields[run]["Amplification"];
    double induction = runs_and_fields[run]["Induction"];
    int tstart = 0;
    int tend = 0;
    
    TChain chain("analysistree/anatree");
    string foldername = path+to_string(run)+"/";
    int n_subruns = glob(string(foldername+"*").data()).size();
    if( n_subruns == 0){
      #if verbose
      cout << foldername << " not found or empty" << endl;
      #endif
      continue;
    }
    #if verbose
    cout << "Processing " << n_subruns << " subruns in run " << run << endl;
    #endif

    string file = "";
    for(int i=0; i< n_subruns; i++){
      
      file=to_string(run)+"-"+to_string(i);
      if(file == "1009-21"){continue;}
      if(path_311data.find("Feb") != string::npos){
        file = file + "-Parser.root";
      }
      else if(path_311data.find("June") != string::npos){
        file = file + "-RecoFast-Parser.root";
      }
      else{
        cout << "ERROR: unknown reconstruction version" << endl;
        return;
      }
      filename=foldername+file;
      if( ExistTest(filename) ){
        #if verbose
        cout << "Adding file: " << filename << endl;
        #endif
        chain.Add(filename.data());
      }
      else{
        #if verbose
        cout << "File " << filename << " not found." << endl;
        #endif
        continue;
      }
    }
    #if verbose
    cout << "Files added. " << endl;
    #endif
    
    int NEntries = (int)chain.GetEntries();
    int  tEventTimeSeconds;
    chain.SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
    chain.GetEntry(0);
    tstart = tEventTimeSeconds;
    chain.GetEntry(NEntries-1);
    tend = tEventTimeSeconds;
    TMyFileHeader *header = new TMyFileHeader(run, drift, amplification, extraction, induction, tstart ,tend);
    TFile ofile(ofilename.data(), "RECREATE");
    if(!ofile.IsOpen()){
      cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
      continue;
    }
    #if verbose
    cout << " Header header_" << to_string(run).data() << " saved in " << ofilename << " has been created " << endl;
    #endif
    
    header->Write(string("header_"+to_string(run)).data());
    ofile.Close();
    
  }//for run
}
