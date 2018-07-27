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


void purity(int run = 840, string cut_type = "before_cuts", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  init_gain_corrections();
  if(gain_corrections.find(run) != gain_corrections.end()){corr = gain_corrections[run];}
  
  gStyle->SetOptFit(1111);
  string path = SelectTrack_Output;
  if(!load_run_lists()){return;}
  //****************************************************************************
  //variables definition here
  float drift = runs_and_fields[run]["Drift"];
  float amplification = runs_and_fields[run]["Extraction"];
  float extraction = runs_and_fields[run]["Amplification"];
  float induction = runs_and_fields[run]["Induction"];
  
//Generated with python macro (all runs)
  
  string filename_nonfitted = dQds_Output+cut_type+"/"+to_string(run)+".root";
  string filename_fitted = dQds_Output+cut_type+"/"+to_string(run)+"_fitted.root";
  string ifilename = "";
  bool read_fit = false;

  if( ExistTest(filename_fitted) ){
    #if verbose
    cout << "      Opening file: " << filename_fitted << endl;
    #endif
    ifilename = filename_fitted;
    read_fit = true;
  }
  else if( ExistTest(filename_nonfitted) ){
    #if verbose
    cout << "      Opening file: " << filename_nonfitted << endl;
    #endif
    ifilename = filename_nonfitted;
  }
  else{
    #if verbose
    cout << "    ERROR: File " << filename_nonfitted << " not found." << endl;
    #endif
    return;
  }
  
  load_cosmics();
  
  
  #if verbose
  cout << "Doing purity scan for run " << run << endl;
  #endif
    
  //Minimum number of tracks to have correct statistics
  int min_number_of_hits = 0;
  //multigraph for mpv, containing view 0 and view 1
  vector<TGraphErrors> purity;
  map<int, vector<TGraphErrors> > purity_ByLEMs;
  
  #if verbose 
  cout << "Initialising graphs..." << endl;
  #endif
  bool are_graph_ok = init_graph_purity_and_charging_up(purity, purity_ByLEMs, "purity");
  if( !are_graph_ok){
    cout << "Error while initialising graph" << endl;
    return;
  }
  
  TFile runfile(ifilename.data(),"READ");

  TH1D* hdQds = 0;
  vector<double> f = {-1,-1};
  
  #if verbose 
  cout << "  Reading histograms..." << endl;
  #endif
  for( int X = 0; X < (int)(100./dx); X++ ){
    int x = dx*X;
    for( auto lem : lems ){
      runfile.GetObject(string("dQds_dx_"+to_string(x)+"_LEM_"+to_string(lem)+"_view0").data(), hdQds);
      if(read_fit){
        f = ReadFit(hdQds,&purity_ByLEMs[lem][0], "", x);
      }
      else{
        f = fit_dQds(hdQds, false, min_number_of_hits, 10000, 1, &purity_ByLEMs[lem][0], "", x);
      }
      runfile.GetObject(string("dQds_dx_"+to_string(x)+"_LEM_"+to_string(lem)+"_view1").data(), hdQds);
      if(read_fit){
        f = ReadFit(hdQds,&purity_ByLEMs[lem][1], "", x);
      }
      else{
        f = fit_dQds(hdQds, false, min_number_of_hits, 10000, 1, &purity_ByLEMs[lem][1], "", x);
      }
    }
    runfile.GetObject(string("dQds_dx_"+to_string(x)+"_view0").data(), hdQds);
    if(read_fit){
      f = ReadFit(hdQds,(&purity[0]), "", x);
    }
    else{
      f = fit_dQds(hdQds, false, min_number_of_hits, 10000, 1, &purity[0], "", x);
    }
    runfile.GetObject(string("dQds_dx_"+to_string(x)+"_view1").data(), hdQds);
    if(read_fit){
      f = ReadFit(hdQds,(&purity[1]), "", x);
    }
    else{
      f = fit_dQds(hdQds, false, min_number_of_hits, 10000, 1, &purity[1], "", x);
    }
  }
  #if verbose
  cout << "    done reading histograms and filling graphs" << endl;
  #endif
  delete hdQds;
  runfile.Close();
  
  string outpath = Purity_Output;
  check_and_mkdir(outpath);
  string outpath = Purity_Output + cut_type;
  check_and_mkdir(outpath);
  string outpath = Purity_Output + "/" + to_string(run) + "/";
  check_and_mkdir(outpath);
  string outfile = outpath + "purity.root"
  TFile ofile(outfile.data(), "RECREATE");
  
  #if verbose
  cout << "    Writing file: " << outfile << endl;
  #endif
  ofile.cd();
  
  vector<double> tau;
  TCanvas can("can","can",1500,750);
  tau.push_back(e_lifetime((&purity[0])));
  tau.push_back(e_lifetime((&purity[1])));
  purity[0].SetTitle(string(string(purity[0].GetTitle())+" \n e lifetime v0: "+to_string_with_precision(tau[0]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[0])+" ppb \n e lifetime v1: "+to_string_with_precision(tau[1]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[1])+" ppb").data());
  purity[1].SetTitle(string(string(purity[1].GetTitle())+" \n e lifetime v0: "+to_string_with_precision(tau[0]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[0])+" ppb \n e lifetime v1: "+to_string_with_precision(tau[1]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[1])+" ppb").data());
  if(TMath::MaxElement(purity[0].GetN(),purity[0].GetY()) > TMath::MaxElement(purity[1].GetN(),purity[1].GetY())){
    purity[0].Draw("AP");
    purity[1].Draw("SAMEP");
  }
  else{
    purity[1].Draw("AP");
    purity[0].Draw("SAMEP");
  }
  ofile.cd();
  purity[0].Write();
  purity[1].Write();
  can.Modified();
  can.Update();
  ofile.cd();
  can.Write();
  
  for( auto lem : lems ){
    #if verbose
    cout << "  Writing purity graph for LEM " << lem << endl;
    #endif
    tau = {};
    TCanvas can_lem(string("can_lem_"+to_string(lem)).data(),string("can_lem_"+to_string(lem)).data(),1500,750);
    tau.push_back(e_lifetime(purity_ByLEMs[lem][0]));
    tau.push_back(e_lifetime(purity_ByLEMs[lem][1]));
    purity_ByLEMs[lem][0].SetTitle(string(string(purity_ByLEMs[lem][0].GetTitle())+" \n e lifetime v0: "+to_string_with_precision(tau[0]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[0])+" ppb \n e lifetime v1: "+to_string_with_precision(tau[1]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[1])+" ppb").data());
    purity_ByLEMs[lem][1].SetTitle(string(string(purity_ByLEMs[lem][1].GetTitle())+" \n e lifetime v0: "+to_string_with_precision(tau[0]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[0])+" ppb \n e lifetime v1: "+to_string_with_precision(tau[1]/1000,3)+" ms. Impurities : "+to_string_with_precision(300/tau[1])+" ppb").data());
    if(TMath::MaxElement(purity_ByLEMs[lem][0].GetN(),purity_ByLEMs[lem][0].GetY()) > TMath::MaxElement(purity_ByLEMs[lem][1].GetN(),purity_ByLEMs[lem][1].GetY())){
      purity_ByLEMs[lem][0].Draw("AP");
      purity_ByLEMs[lem][1].Draw("SAMEP");
    }
    else{
      purity_ByLEMs[lem][1].Draw("AP");
      purity_ByLEMs[lem][0].Draw("SAMEP");
    }
    ofile.cd();
    purity_ByLEMs[lem][0].Write();
    purity_ByLEMs[lem][1].Write();
    can_lem.Modified();
    can_lem.Update();
    ofile.cd();
    can_lem.Write();
  }
  #if verbose
  cout << "all done!" << endl;
  #endif
  ofile.Close();
}//end macro
