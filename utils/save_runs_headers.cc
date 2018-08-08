#include <stdlib.h>
#include <stdio.h>
#include <311Lib.h>
#include <TMyFileHeader.h>

using namespace std;

void save_runs_headers(vector<int> run_list = {}, string version = "Feb"){
  if(!Load_Version(version)){return false;}
  if(!load_run_lists()){return false;}
  
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
  int failed = 0;
  
  for(auto run : run_list){
  
    if(!load_rho_run(run)){failed++;continue;}
    string ofilename = path_output+to_string(run)+".root";
    double drift = runs_and_fields[run]["Drift"];
    double amplification = runs_and_fields[run]["Amplification"];
    double extraction = runs_and_fields[run]["Extraction"];
    double induction = runs_and_fields[run]["Induction"];
    double extr_eff_simu = get_eff(extraction, amplification);
    double extr_eff_gushin = get_eff(extraction);
    double ind_eff = get_eff(-1,amplification,induction);
    double density_correction_factor = gain_correction_for_rho(rho_run, amplification);
    double drift_correction_factor = correct_for_drift(drift);
    int tstart = 0;
    int tend = 0;
    
    TChain chain("analysistree/anatree");
    string foldername = path+to_string(run)+"/";
    int n_subruns = glob(string(foldername+"*").data()).size();
    if( n_subruns == 0){
      #if verbose
      cout << foldername << " not found or empty" << endl;
      #endif
      failed++;
      continue;
    }
    #if verbose
    cout << "Processing " << n_subruns << " subruns in run " << run << endl;
    #endif

    string file = "";
    for(int i=0; i< n_subruns; i++){
      
      file=to_string(run)+"-"+to_string(i);
      if(file == "1009-21"){failed++;continue;}
      if(path_311data.find("Feb") != string::npos){
        file = file + "-Parser.root";
      }
      else if(path_311data.find("June") != string::npos){
        file = file + "-RecoFast-Parser.root";
      }
      else{
        cout << "ERROR: unknown reconstruction version" << endl;
        return false;
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
        failed++;
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
    string command  = path_311analysis+"slow_control/Pcotte_getdata.sh " + to_string(tstart) + " " + to_string(tend) + " " + to_string(run);
    system(command.data());
    if(!ExistTest(slow_control + "/" + to_string(run) + "/" + to_string(run) + ".txt")){
      cout << "ERROR: can not get slow control data for run " << run << endl;
      failed++;
      continue;
    }
    vector<double> rho = pressure(to_string(run));
    if(rho[0] < 0){
      cout << "ERROR: can not get pressure and temperature for run " << run << endl;
      failed++;
      continue;
    }
    
    TMyFileHeader *header = new TMyFileHeader(run, drift, amplification, extraction, induction, tstart ,tend,0, rho[0], rho[1], extr_eff_simu, extr_eff_gushin, ind_eff, density_correction_factor, drift_correction_factor);
    TFile ofile(ofilename.data(), "RECREATE");
    if(!ofile.IsOpen()){
      cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
      failed++;
      continue;
    }
    #if verbose
    cout << " Header header_" << to_string(run).data() << " saved in " << ofilename << " has been created " << endl;
    #endif
    
    header->Write(string("header_"+to_string(run)).data());
    ofile.Close();
    
  }//for run
  if(failed == run_list.size()){
    cout << "ERROR in save_run_headers: all runs failed" << endl;
    return false;
  }
  return true;
}
