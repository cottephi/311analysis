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
#include <stdlib.h>
#include <stdio.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}


void hits_goodness_of_fit_vs_dqds(int run = 840, string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_dQ = m_dQ;
  method_ds = m_ds;
  if(!Load_Version(version)){return;}

  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  
  #if verbose 
  cout << "Initialising histos..." << endl;
  #endif
  vector<TH2D> HGoF;
  HGoF.push_back(TH2D("HGoF_0","HGoF_0",100,0,50,100,-2,10));
  HGoF.push_back(TH2D("HGoF_1","HGoF_1",100,0,50,100,-2,10));
  HGoF.push_back(TH2D("HGoF_summed_views","HGoF_summed_views",100,0,50,100,-2,10));
  map<int, vector<TH2D>  > HGoF_ByLEMs;
  vector<TH2D> HGoF_integral_local;
  HGoF_integral_local.push_back(TH2D("HGoF_integral_local_integral_local_0","HGoF_integral_local_0",100,0,50,100,-2,10));
  HGoF_integral_local.push_back(TH2D("HGoF_integral_local_1","HGoF_integral_local_1",100,0,50,100,-2,10));
  HGoF_integral_local.push_back(TH2D("HGoF_integral_local_summed_views","HGoF_integral_local_summed_views",100,0,50,100,-2,10));
  map<int, vector<TH2D>  > HGoF_integral_local_ByLEMs;
  for(auto lem : lems){
    HGoF_ByLEMs[lem].push_back(TH2D(string("HGoF_0_LEM_"+to_string(lem)).data(),string("HGoF_0_LEM_"+to_string(lem)).data(),100,0,50,100,-2,10));
    HGoF_ByLEMs[lem].push_back(TH2D(string("HGoF_1_LEM_"+to_string(lem)).data(),string("HGoF_1_LEM_"+to_string(lem)).data(),100,0,50,100,-2,10));
    HGoF_ByLEMs[lem].push_back(TH2D(string("HGoF_summed_views_LEM_"+to_string(lem)).data(),string("HGoF_summed_views_LEM_"+to_string(lem)).data(),100,0,50,100,-2,10));
    
    HGoF_integral_local_ByLEMs[lem].push_back(TH2D(string("HGoF_integral_local_0_LEM_"+to_string(lem)).data(),string("HGoF_integral_local_0_LEM_"+to_string(lem)).data(),100,0,50,100,-2,10));
    HGoF_integral_local_ByLEMs[lem].push_back(TH2D(string("HGoF_integral_local_1_LEM_"+to_string(lem)).data(),string("HGoF_integral_local_1_LEM_"+to_string(lem)).data(),100,0,50,100,-2,10));
    HGoF_integral_local_ByLEMs[lem].push_back(TH2D(string("HGoF_integral_local_summed_views_LEM_"+to_string(lem)).data(),string("HGoF_integral_local_summed_views_LEM_"+to_string(lem)).data(),100,0,50,100,-2,10));
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
      HGoF[h.view].Fill(h.dq_sum/h.dx_3D,h.GoF);
      HGoF_integral_local[h.view].Fill(h.dq_integral/h.dx_local,h.GoF);
      HGoF_ByLEMs[h.lem][h.view].Fill(h.dq_sum/h.dx_3D,h.GoF);
      HGoF_integral_local_ByLEMs[h.lem][h.view].Fill(h.dq_integral/h.dx_local,h.GoF);
    }
  }

  #if verbose
  cout << "Plotting hits goodness of fit vs ds for run " << run << endl;
  #endif
  string outpath = path_311data+"divers/";
  check_and_mkdir(outpath);
  outpath = outpath + to_string(run)+"/";
  check_and_mkdir(outpath);
  TFile ofile(string(outpath+"HGoF_vs_dqds.root").data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << " ERROR: TFile " << outpath << "HGoF_vs_dqds.root can't be created " << endl;
    return;
  }
  #if verbose
  cout << " TFile " << outpath << "HGoF_vs_dqds.root has been created " << endl;
  #endif
  
  for(auto his : HGoF){
    his.Write();
    his.SetMaximum(5000);
    his.SetMinimum(10);
    his.GetXaxis()->SetTitle("dQds");
    his.GetYaxis()->SetTitle("GoF");
    his.Draw("COLZ");
    gPad->SetName(string(outpath+string(his.GetName())+"_pad").data());
    gPad->SaveAs(string(outpath+string(his.GetName())+".png").data());
    gPad->Write();
  }
  for(auto his : HGoF_integral_local){
    his.Write();
    his.SetMaximum(5000);
    his.SetMinimum(10);
    his.GetXaxis()->SetTitle("dQds");
    his.GetYaxis()->SetTitle("GoF");
    his.Draw("COLZ");
    gPad->SetName(string(outpath+string(his.GetName())+"_pad").data());
    gPad->SaveAs(string(outpath+string(his.GetName())+".png").data());
    gPad->Write();
  }
  for( auto lem : lems ){
    for(auto his : HGoF_ByLEMs[lem]){
      his.Write();
      his.SetMaximum(1000);
      his.SetMinimum(10);
      his.GetXaxis()->SetTitle("dQds");
      his.GetYaxis()->SetTitle("GoF");
      his.Draw("COLZ");
      gPad->SetName(string(outpath+string(his.GetName())+"_pad").data());
      gPad->SaveAs(string(outpath+string(his.GetName())+".png").data());
      gPad->Write();
    }
    for(auto his : HGoF_integral_local_ByLEMs[lem]){
      his.Write();
      his.SetMaximum(1000);
      his.SetMinimum(10);
      his.GetXaxis()->SetTitle("dQds");
      his.GetYaxis()->SetTitle("GoF");
      his.Draw("COLZ");
      gPad->SetName(string(outpath+string(his.GetName())+"_pad").data());
      gPad->SaveAs(string(outpath+string(his.GetName())+".png").data());
      gPad->Write();
    }
  }
  ofile.Close();
}//end macro
