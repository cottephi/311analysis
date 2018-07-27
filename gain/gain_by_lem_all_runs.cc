////////////////////////////////////////////////////////////////////////////////
//
// Produce gain curve using the MPV of the dQdx distribution
// Takes as input the root file produced by utils/track_cuts.cc
// Produce a root file with all the dQdx distributions for every Electric fields
// and a TGraphErrors with the gain curve
//
//TODO: substitute to_string() in histname
//TODO: add scan type selection (if induction, extraction, etc...)
//
// Usage:
// Toggle the scan_type (Extraction, Induction, Amplification) and scan_num (1,2)
// to run the macro and producte the MPV vs field of the scan
//
//
////////////////////////////////////////////////////////////////////////////////

#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <TMyFileHeader.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}

void gain_by_lem_all_runs(vector<int> run_list = {}, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  gErrorIgnoreLevel = kError;
  
  string path = gain_by_lem_Output + cut_type + "/";
  if(run_list.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path << "*..." << endl;
    #endif
    string wildcard_path = path + "*";
    for( auto irun : glob(wildcard_path) ){
      if( string(irun).find("std_dev") != string::npos ){continue;}
      run_list.push_back(atoi(irun.data()));
    }
  }
  
//  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);
  
  TH1D stddev("std_dev","std_dev",20,0,0.4);
  TH1D stddev_view0("std_dev_view0","std_dev_view0",20,0,0.4);
  TH1D stddev_view1("std_dev_view1","std_dev_view1",20,0,0.4);

  load_cosmics();
  
  for(auto run : run_list){
    #if verbose
    cout << "  Processing run " << run << endl;
    #endif
    string ifilename = path + to_string(run) + "/gain.root";
    TFile ifile(ifilename.data(),"READ");
    TH1D *distri = 0;
    ifile.GetObject("distri_ByLEMs",distri);
    stddev.Fill(distri->GetStdDev());
    distri = 0;
    ifile.GetObject("distri_ByLEMs_view0",distri);
    stddev_view0.Fill(distri->GetStdDev());
    distri = 0;
    ifile.GetObject("distri_ByLEMs_view1",distri);
    stddev_view1.Fill(distri->GetStdDev());
    delete distri; distri = 0;
    ifile.Close();
  }
  
  
  TFile ofile(string(path+"std_dev_all_run.root").data(),"RECREATE");
  #if verbose
  cout << "Writing file " << ofile.GetName() << endl;
  #endif
  
  ofile.cd();
  stddev.Draw();
  gPad->SetName(string(string(stddev.GetName()) + "_pad").data());
  gPad->Write();
  gPad->SaveAs(string(path + string(stddev.GetName()) + ".png").data());
  
  ofile.cd();
  stddev_view0.Draw();
  gPad->SetName(string(string(stddev_view0.GetName()) + "_pad").data());
  gPad->Write();
  gPad->SaveAs(string(path + string(stddev_view0.GetName()) + ".png").data());
  
  ofile.cd();
  stddev_view1.Draw();
  gPad->SetName(string(string(stddev_view1.GetName()) + "_pad").data());
  gPad->Write();
  gPad->SaveAs(string(path + string(stddev_view1.GetName()) + ".png").data());
  
  delete gPad;
  ofile.Close();
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
