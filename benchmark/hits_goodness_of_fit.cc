////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <cstdlib>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}


void hits_goodness_of_fit(int run = 840, string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_dQ = m_dQ;
  method_ds = m_ds;
  if(!Load_Version(version)){return;}

  gErrorIgnoreLevel = kError;
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.35);                
  // Set y-position (fraction of pad size)
  gStyle->SetStatX(0.85);                
  // Set x-position (fraction of pad size)
  gStyle->SetStatW(0.2);                
  // Set width of stat-box (fraction of pad size)
  gStyle->SetStatH(0.15);                
  // Set height of stat-box (fraction of pad size)
  string path = SelectTrack_Output;
  if(!load_run_lists()){return;}
  vector<string> run_list;
  load_cosmics();
  int min_number_of_hits = 0;
  int time = 1;

  //****************************************************************************
  //variables definition here
  
//Generated with python macro (all runs)
  
  
  
  //multigraph for mpv, containing view 0 and view 1
  vector<TH1D> HGoF;
  HGoF.push_back(TH1D("HGoF_0","HGoF_0",100,0,10));
  HGoF.push_back(TH1D("HGoF_1","HGoF_1",100,0,10));
  HGoF.push_back(TH1D("HGoF_summed_views","HGoF_summed_views",100,0,10));
  map<int, vector<TH1D>  > HGoF_ByLEMs;
  #if verbose 
  cout << "Initialising histos..." << endl;
  #endif
  for(auto lem : lems){
    HGoF_ByLEMs[lem].push_back(TH1D(string("HGoF_0_LEM_"+to_string(lem)).data(),string("HGoF_0_LEM_"+to_string(lem)).data(),100,0,10));
    HGoF_ByLEMs[lem].push_back(TH1D(string("HGoF_1_LEM_"+to_string(lem)).data(),string("HGoF_1_LEM_"+to_string(lem)).data(),100,0,10));
    HGoF_ByLEMs[lem].push_back(TH1D(string("HGoF_summed_views_LEM_"+to_string(lem)).data(),string("HGoF_summed_views_LEM_"+to_string(lem)).data(),100,0,10));
  }
  #if verbose 
  cout << "...done" << endl;
  #endif
  
  string inputpath = SelectTrack_Output + "/before_cuts/tracks_" + to_string(run) + ".root";
  TFile ifile(inputpath.data(), "READ");
  if(!ifile.IsOpen()){
    cout << "ERROR: can not open " << inputpath << endl;
    return;
  }
  vector<track> tracks;
  track *single_track = 0;
  TTree *tree = 0;
  ifile.GetObject("tracks", tree);
  tree->SetBranchAddress("track", &single_track);
  #if verbose 
  cout << "Loading tracks..." << endl;
  #endif
  for(int i = 0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    tracks.push_back(*single_track);
  }
  #if verbose 
  cout << "...done" << endl;
  #endif
  for(auto tr : tracks){
    for(auto h : tr.hits_trk){
      if(!isGood_lem(h.lem)){continue;}
      if(!IsGoodChan(h)){continue;}
      HGoF[h.view].Fill(h.GoF);
      HGoF_ByLEMs[h.lem][h.view].Fill(h.GoF);
    }
  }
  
  #if verbose
  cout << "Plotting hits goodness of fit for run " << run << endl;
  #endif
  string outpath = path_311data+"divers/";
  check_and_mkdir(outpath);
  outpath = outpath + to_string(run)+"/";
  check_and_mkdir(outpath);
  TFile ofile(string(outpath+"HGoF.root").data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << " ERROR: TFile " << outpath << "HGoF.root can't be created " << endl;
    return;
  }
  #if verbose
  cout << " TFile " << outpath << "HGoF.root has been created " << endl;
  #endif
  
  HGoF[0].Write();
  HGoF[1].Write();
  for( auto lem : lems ){
    HGoF_ByLEMs[lem][0].Write();
    HGoF_ByLEMs[lem][1].Write();
  }
  ofile.Close();
}//end macro
