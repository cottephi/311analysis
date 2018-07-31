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


void charging_up(int run = 840, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  gErrorIgnoreLevel = kError;
  
  string corrected_or_not = "";
  if(cut_type.find("tg") != string::npos){corrected_or_not = "_Dx_Corrected";}
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
  
  string inpath = dQds_charging_up_Output + cut_type + "/" + to_string(run) + "/";
  for( auto irun : glob(inpath+"*fitted*") ){
    run_list.push_back(irun);
  }
  int Nfitted = run_list.size();
  run_list.clear();
  for( auto irun : glob(inpath+"*") ){
    run_list.push_back(irun);
  }
  int Nfiles = run_list.size();
  int Nsubruns = Nfiles - Nfitted;
  
  if(Nsubruns == 0){
    cout << "Error: no subruns found in " << inpath << endl;
    cout << "ls in this folder: " << endl;
    system(string("ls -l " + inpath).data());
    return;
  }
  
  #if verbose
  cout << "Doing charging up analysis for run " << run << endl;
  #endif
  
  
  //multigraph for mpv, containing view 0 and view 1
  vector<TGraphErrors> mg_charging_up;
  map<int, vector<TGraphErrors>  > mg_charging_up_ByLEMs;
  
  #if verbose 
  cout << "Initialising graphs..." << endl;
  #endif
  bool are_graph_ok = init_graph_purity_and_charging_up(mg_charging_up, mg_charging_up_ByLEMs);
  if( !are_graph_ok){
    cout << "Error while initialising graph" << endl;
    return;
  }
  bool found_first = false;
  bool recreate_fit_file = true;
  double initgain0 = 0;
  double initgain1 = 0;
  int time0 = 0;
  
  for(int subrun = 0; subrun < Nsubruns; subrun++){
  
    string filename_nonfitted = inpath+to_string(subrun)+".root";
    string filename_fitted = inpath+to_string(subrun)+"_fitted.root";
    int files_not_found = 0;

    TFile *runfile_fitted = new TFile();
    bool read_fit = false;
    string ifilename = "";
    int i_read_fit = read_or_do_fit(filename_fitted, filename_nonfitted, recreate_fit_file, ifilename, files_not_found);
    if(i_read_fit == 3){return;}
    if(i_read_fit == 1){read_fit = true;}
    TFile runfile(ifilename.data(),"READ");
    if( !read_fit ){
      #if verbose
      cout << "    Recording fitted histograms in " << filename_fitted << endl;
      #endif
      runfile_fitted = TFile::Open(filename_fitted.data(), "RECREATE");
      if(!runfile_fitted->IsOpen()){
        cout << "ERROR: could not open " << filename_fitted << " for writing." << endl;
        return;
      }
    }
  
    TH1D* hdQds = 0;
    vector<double> f0 = {-1,-1};
    vector<double> f1 = {-1,-1};
    
    #if verbose 
    cout << "  Reading dQds histograms of run " << run << "..." << endl;
    #endif
    runfile.GetObject(string("dQds_view0" + corrected_or_not).data(), hdQds);
    time = atoi(hdQds->GetTitle());
    if(read_fit){
      f0 = ReadFit(hdQds);
    }
    else{
      f0 = fit_dQds(hdQds, false, min_number_of_hits, 10000, 0.25);
      runfile_fitted->cd();
      hdQds->Write();
    }
    runfile.GetObject(string("dQds_view1" + corrected_or_not).data(), hdQds);
    if(read_fit){
      f1 = ReadFit(hdQds);
    }
    else{
      f1 = fit_dQds(hdQds, false, min_number_of_hits, 10000, 0.25);
      runfile_fitted->cd();
      hdQds->Write();
    }
    if(mpv_cosmics > 0){
      f0[0] = f0[0]/mpv_cosmics;
      f0[1] = f0[1]/mpv_cosmics;
      f1[0] = f1[0]/mpv_cosmics;
      f1[1] = f1[1]/mpv_cosmics;
    }
    if(f1[0] > 0 && f0[0] > 0){
      if(found_first == false){
        time0 = time;
        initgain0 = f0[0];
        initgain1 = f1[0];
        mg_charging_up[0].SetPoint(mg_charging_up[0].GetN(),time-time0,f0[0]);
        mg_charging_up[0].SetPointError(mg_charging_up[0].GetN()-1,0,f0[1]);
        mg_charging_up[1].SetPoint(mg_charging_up[1].GetN(),time-time0,f1[0]);
        mg_charging_up[1].SetPointError(mg_charging_up[1].GetN()-1,0,f1[1]);
        mg_charging_up[2].SetPoint(mg_charging_up[2].GetN(),time-time0,f0[0]+f1[0]);
        mg_charging_up[2].SetPointError(mg_charging_up[2].GetN()-1,0,TMath::Sqrt(pow(f0[1],2)+pow(f1[1],2)));
        found_first = true;
      }
      else if(initgain0 > 0 && initgain1 > 0){
        if(f0[0] > 0.7*initgain0 && f1[0] > 0.7*initgain1){
          mg_charging_up[0].SetPoint(mg_charging_up[0].GetN(),time-time0,f0[0]);
          mg_charging_up[0].SetPointError(mg_charging_up[0].GetN()-1,0,f0[1]);
          mg_charging_up[1].SetPoint(mg_charging_up[1].GetN(),time-time0,f1[0]);
          mg_charging_up[1].SetPointError(mg_charging_up[1].GetN()-1,0,f1[1]);
          mg_charging_up[2].SetPoint(mg_charging_up[2].GetN(),time-time0,f0[0]+f1[0]);
          mg_charging_up[2].SetPointError(mg_charging_up[2].GetN()-1,0,TMath::Sqrt(pow(f0[1],2)+pow(f1[1],2)));
        }
      }
    }
    for( auto lem : lems ){
      runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view0" + corrected_or_not).data(), hdQds);
      if(read_fit){
        f0 = ReadFit(hdQds);
      }
      else{
        f0 = fit_dQds(hdQds, false, min_number_of_hits, 10000, 1);
        runfile_fitted->cd();
        hdQds->Write();
      }
      runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view1" + corrected_or_not).data(), hdQds);
      if(read_fit){
        f1 = ReadFit(hdQds);
      }
      else{
        f1 = fit_dQds(hdQds, false, min_number_of_hits, 10000, 1);
        runfile_fitted->cd();
        hdQds->Write();
      }
      if(mpv_cosmics > 0){
        f0[0] = f0[0]/mpv_cosmics;
        f0[1] = f0[1]/mpv_cosmics;
        f1[0] = f1[0]/mpv_cosmics;
        f1[1] = f1[1]/mpv_cosmics;
      }
      if(f1[0] > 0 && f0[0] > 0){
        if(found_first == false){
          time0 = time;
          initgain0 = f0[0];
          initgain1 = f1[0];
          mg_charging_up_ByLEMs[lem][0].SetPoint(mg_charging_up_ByLEMs[lem][0].GetN(),time-time0,f0[0]);
          mg_charging_up_ByLEMs[lem][0].SetPointError(mg_charging_up_ByLEMs[lem][0].GetN()-1,0,f0[1]);
          mg_charging_up_ByLEMs[lem][1].SetPoint(mg_charging_up_ByLEMs[lem][1].GetN(),time-time0,f1[0]);
          mg_charging_up_ByLEMs[lem][1].SetPointError(mg_charging_up_ByLEMs[lem][1].GetN()-1,0,f1[1]);
          mg_charging_up_ByLEMs[lem][2].SetPoint(mg_charging_up_ByLEMs[lem][2].GetN(),time-time0,f0[0]+f1[0]);
          mg_charging_up_ByLEMs[lem][2].SetPointError(mg_charging_up_ByLEMs[lem][2].GetN()-1,0,TMath::Sqrt(pow(f0[1],2)+pow(f1[1],2)));
          found_first = true;
        }
        else if(initgain0 > 0 && initgain1 > 0){
          if(f0[0] > 0.7*initgain0 && f1[0] > 0.7*initgain1){
            mg_charging_up_ByLEMs[lem][0].SetPoint(mg_charging_up_ByLEMs[lem][0].GetN(),time-time0,f0[0]);
            mg_charging_up_ByLEMs[lem][0].SetPointError(mg_charging_up_ByLEMs[lem][0].GetN()-1,0,f0[1]);
            mg_charging_up_ByLEMs[lem][1].SetPoint(mg_charging_up_ByLEMs[lem][1].GetN(),time-time0,f1[0]);
            mg_charging_up_ByLEMs[lem][1].SetPointError(mg_charging_up_ByLEMs[lem][1].GetN()-1,0,f1[1]);
            mg_charging_up_ByLEMs[lem][2].SetPoint(mg_charging_up_ByLEMs[lem][2].GetN(),time-time0,f0[0]+f1[0]);
            mg_charging_up_ByLEMs[lem][2].SetPointError(mg_charging_up_ByLEMs[lem][2].GetN()-1,0,TMath::Sqrt(pow(f0[1],2)+pow(f1[1],2)));
          }
        }
      }
    }
    #if verbose
    cout << "    done reading histograms and filling graphs" << endl;
    #endif
    delete hdQds;
    if(runfile_fitted->IsOpen()){
      runfile_fitted->Close();
    }
    delete runfile_fitted;
    runfile.Close();
  }//subruns
  
  
//attempt to fit if histos are not emtpy
      
  sum_views(mg_charging_up, mg_charging_up_ByLEMs);

  string outpath = charging_up_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type + "/";
  check_and_mkdir(outpath);
  outpath = outpath + "/"+to_string(run) + "/";
  check_and_mkdir(outpath);
  string outfile = outpath+"charging_up.root";
  TFile ofile(outfile.data(), "RECREATE");
  #if verbose
  cout << "Writing file: " << string(ofile.GetName()).data() << endl;
  #endif
  fit_charging_up(mg_charging_up[2]);
  fit_charging_up(mg_charging_up[1]);
  fit_charging_up(mg_charging_up[0]);
  mg_charging_up[2].Draw("AP");
  mg_charging_up[2].SetMinimum(0);
  mg_charging_up[2].GetXaxis()->SetTimeDisplay(1);
  mg_charging_up[2].GetXaxis()->SetTimeFormat("%Hh%M");
  mg_charging_up[2].Draw("AP"); mg_charging_up[0].Draw("SAMEP"); mg_charging_up[1].Draw("SAMEP");
  ofile.cd();
  mg_charging_up[0].Write(); mg_charging_up[1].Write(); mg_charging_up[2].Write();
  gPad->Modified();
  gPad->Update();
  gPad->SetName("charging_up_pad");
  ofile.cd();
  gPad->SaveAs(string(outpath + string(gPad->GetName())+".png").data());
  gPad->Write();
  delete gPad;
  
  for( auto lem : lems ){
    #if verbose
    cout << "  Writing charging up graph for LEM " << lem << endl;
    #endif
    fit_charging_up(mg_charging_up_ByLEMs[lem][2]);
    fit_charging_up(mg_charging_up_ByLEMs[lem][1]);
    fit_charging_up(mg_charging_up_ByLEMs[lem][0]);
    mg_charging_up_ByLEMs[lem][2].Draw("AP");
    mg_charging_up_ByLEMs[lem][2].SetMinimum(0);
    mg_charging_up_ByLEMs[lem][2].GetXaxis()->SetTimeDisplay(1);
    mg_charging_up_ByLEMs[lem][2].GetXaxis()->SetTimeFormat("%Hh%M");
    mg_charging_up_ByLEMs[lem][2].Draw("AP"); mg_charging_up_ByLEMs[lem][0].Draw("SAMEP"); mg_charging_up_ByLEMs[lem][1].Draw("SAMEP");
    ofile.cd();
    mg_charging_up_ByLEMs[lem][0].Write();mg_charging_up_ByLEMs[lem][1].Write();mg_charging_up_ByLEMs[lem][2].Write();
    gPad->Modified();
    gPad->Update();
    gPad->SetName(string("chargin_up_LEM_"+ to_string(lem)+"_pad").data());
    ofile.cd();
    gPad->Write();
    gPad->SaveAs(string(outpath + string(gPad->GetName())+".png").data());
    delete gPad;
  }//lems
  #if verbose
  cout << "all done!" << endl;
  #endif
  ofile.Close();
}//end macro
