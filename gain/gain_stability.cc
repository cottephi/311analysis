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


void gain_stability(const char* c_cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  gErrorIgnoreLevel = kWarning;
//  gStyle->SetOptStat(0);
  
  bool recreate_fit_file = true;
  double bin_width = 0.05;
  double gain_max = 4.;
  int nbins = gain_max/bin_width;
  double max = 0;
  
  string cut_type = c_cut_type;

  if(!load_run_lists()){return;}
  
  //****************************************************************************
  //variables definition here
  
//Generated with python macro (all runs), 50% drift tolerence
  vector<int> scan_nums = {/*16*/};
  for(int i = 1; i < 27; i++){
    scan_nums.push_back(i);
  }
  //system(string("rm "+GainAna_Output+"*.root").data());
  load_cosmics();
  
  string scan_type = "Identical";

  int empty_scan_num = 0;
  for( auto scan_num : scan_nums ){

    string histname = "mpv_run_"+to_string(scan_num);
    TH1D* mpv_run = new TH1D(histname.data(),histname.data(), nbins, 0, gain_max);
    histname = "mpv_run_"+to_string(scan_num)+"_view0";
    TH1D* mpv_run_view0 = new TH1D(histname.data(),histname.data(), nbins, 0, gain_max);
    histname = "mpv_run_"+to_string(scan_num)"_view1";
    TH1D* mpv_run_view1 = new TH1D(histname.data(),histname.data(), nbins, 0, gain_max);
    map<int, TH1D*> mpv_run_ByLEMs;
    map<int, TH1D*> mpv_run_ByLEMs_view0;
    map<int, TH1D*> mpv_run_ByLEMs_view1;
    histname = "mpv_run_"+to_string(scan_num)+"_Dx_Corrected";
    TH1D* mpv_run_Dx_Corrected = new TH1D(histname.data(),histname.data(), nbins, 0, gain_max);
    histname = "mpv_run_"+to_string(scan_num)+"_Dx_Corrected_view0";
    TH1D* mpv_run_Dx_Corrected_view0 = new TH1D(histname.data(),histname.data(), nbins, 0, gain_max);
    histname = "mpv_run_"+to_string(scan_num)+"_Dx_Corrected_view1";
    TH1D* mpv_run_Dx_Corrected_view1 = new TH1D(histname.data(),histname.data(), nbins, 0, gain_max);
    map<int, TH1D*> mpv_run_ByLEMs_Dx_Corrected;
    map<int, TH1D*> mpv_run_ByLEMs_Dx_Corrected_view0;
    map<int, TH1D*> mpv_run_ByLEMs_Dx_Corrected_view1;
    #if verbose 
    cout << "  Initialising stability histos..." << endl;
    #endif

    for( auto lem : lems ){
      histname = "mpv_run_LEM_"+to_string(lem)+"_"+to_string(scan_num);
      mpv_run_ByLEMs[lem] = new TH1D(histname.data(), histname.data(), nbins, 0, gain_max);
      histname = "mpv_run_LEM_"+to_string(lem)+"_"+to_string(scan_num)+"_view0";
      mpv_run_ByLEMs_view0[lem] = new TH1D(histname.data(), histname.data(), nbins, 0, gain_max);
      histname = "mpv_run_LEM_"+to_string(lem)+"_"+to_string(scan_num)+"_view1";
      mpv_run_ByLEMs_view1[lem] = new TH1D(histname.data(), histname.data(), nbins, 0, gain_max);
      histname = "mpv_run_LEM_"+to_string(lem)+"_"+to_string(scan_num)+"_Dx_Corrected";
      mpv_run_ByLEMs_Dx_Corrected[lem] = new TH1D(histname.data(), histname.data(), nbins, 0, gain_max);
      histname = "mpv_run_LEM_"+to_string(lem)+"_"+to_string(scan_num)+"_Dx_Corrected_view0";
      mpv_run_ByLEMs_Dx_Corrected_view0[lem] = new TH1D(histname.data(), histname.data(), nbins, 0, gain_max);
      histname = "mpv_run_LEM_"+to_string(lem)+"_"+to_string(scan_num)+"_Dx_Corrected_view1";
      mpv_run_ByLEMs_Dx_Corrected_view1[lem] = new TH1D(histname.data(), histname.data(), nbins, 0, gain_max);
    }//for lem
    

    //Minimum number of hits to have correct statistics
    int min_number_of_hits = 200;

    vector<int> run_list;
    vector<float> dummy1;
    vector<string> dummy2;
    vector<float> dummy3;

    if(!load_runs(scan_type, scan_num, run_list, dummy1, dummy2, dummy3)){continue;}
    
    #if verbose
    cout << "  Doing " << scan_type << " scan number " << scan_num << endl;
    #endif
    
    const int n_runs = run_list.size();
    
    int files_not_found = 0;
    for(auto run : run_list){

      string filename_nonfitted = dQds_Output+cut_type+"/"+to_string(run)+".root";
      string filename_fitted = dQds_Output+cut_type+"/"+to_string(run)+"_fitted.root";
      string ifilename = "";
      bool read_fit = false;
      if( ExistTest(filename_fitted) && !recreate_fit_file){
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
        cout << "    File " << filename_nonfitted << " not found." << endl;
        #endif
        files_not_found++;
        continue;
      }
      TFile runfile(ifilename.data(),"READ");
      TFile *runfile_fitted = new TFile();
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
      TString FunName = "";
      int nhits_tot = 0;
      vector<double> f0 = {-1,-1};
      vector<double> f1 = {-1,-1};
      int min_number_of_hits = 200;
        
      #if verbose 
      cout << "  Reading histograms..." << endl;
      #endif
      runfile.GetObject("dQds_view0", hdQds);
      nhits_tot += hdQds->GetEntries();
      if(read_fit){
        f0 = ReadFit(hdQds);
      }
      else{
        f0 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
        runfile_fitted->cd();
        hdQds->Write();
      }
      if(mpv_cosmics > 0 && f0[0]){
        mpv_run_view0->Fill(f0[0]/mpv_cosmics);
        if(f0[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        if(f0[0]/mpv_cosmics > max){max = f0[0]/mpv_cosmics;}
      }
      else if(f0[0] > 0){mpv_run_view0->Fill(f0[0]);}
      runfile.GetObject("dQds_view1", hdQds);
      nhits_tot += hdQds->GetEntries();
      if(read_fit){
        f1 = ReadFit(hdQds);
      }
      else{
        f1 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
        runfile_fitted->cd();
        hdQds->Write();
      }
      if(mpv_cosmics > 0 && f1[0]){
        mpv_run_view1->Fill(f1[0]/mpv_cosmics);
        if(f1[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        if(f1[0]/mpv_cosmics > max){max = f1[0]/mpv_cosmics;}
      }
      else if(f1[0] > 0){mpv_run_view1->Fill(f1[0]);}
      if(mpv_cosmics > 0 && f0[0] > 0 && f1[0] > 0){
        mpv_run->Fill((f0[0]+f1[0])/mpv_cosmics);
        if((f0[0]+f1[0])/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        if((f0[0]+f1[0])/mpv_cosmics > max){max = (f0[0]+f1[0])/mpv_cosmics;}
      }
      else if(f0[0] > 0 && f1[0] > 0){mpv_run->Fill(f0[0]+f1[0]);}
      if( nhits_tot < min_number_of_hits ){
        #if verbose
        cout << "    Insuficient number of hits (" << nhits_tot << ") for run " << run << ". Skipping this run. \nnext run" << endl;
        #endif
        if(runfile_fitted->IsOpen()){
          runfile_fitted->Close();
        }
        delete runfile_fitted;
        runfile.Close();
        continue;
      }
      else{
        #if verbose
        cout << "    Processing " << nhits_tot << " hits for run " << run << endl;
        #endif
      }
      
      runfile.GetObject("dQds_view0_Dx_Corrected", hdQds);
      nhits_tot += hdQds->GetEntries();
      if(read_fit){
        f0 = ReadFit(hdQds);
      }
      else{
        f0 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
        runfile_fitted->cd();
        hdQds->Write();
      }
      if(mpv_cosmics > 0 && f0[0]){
        mpv_run_Dx_Corrected_view0->Fill(f0[0]/mpv_cosmics);
        if(f0[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        if(f0[0]/mpv_cosmics > max){max = f0[0]/mpv_cosmics;}
      }
      else if(f0[0] > 0){mpv_run_Dx_Corrected_view0->Fill(f0[0]);}
      runfile.GetObject("dQds_view1_Dx_Corrected", hdQds);
      nhits_tot += hdQds->GetEntries();
      if(read_fit){
        f1 = ReadFit(hdQds);
      }
      else{
        f1 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
        runfile_fitted->cd();
        hdQds->Write();
      }
      if(mpv_cosmics > 0 && f1[0]){
        mpv_run_Dx_Corrected_view1->Fill(f1[0]/mpv_cosmics);
        if(f1[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        if(f1[0]/mpv_cosmics > max){max = f1[0]/mpv_cosmics;}
      }
      else if(f1[0] > 0){mpv_run_Dx_Corrected_view1->Fill(f1[0]);}
      if(mpv_cosmics > 0 && f0[0] > 0 && f1[0] > 0){
        mpv_run_Dx_Corrected->Fill((f0[0]+f1[0])/mpv_cosmics);
        if((f0[0]+f1[0])/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        if((f0[0]+f1[0])/mpv_cosmics > max){max = (f0[0]+f1[0])/mpv_cosmics;}
      }
      else if(f0[0] > 0 && f1[0] > 0){mpv_run_Dx_Corrected->Fill(f0[0]+f1[0]);}
      
      for( auto lem : lems ){
        runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view0").data(), hdQds);
        if(read_fit){
          f0 = ReadFit(hdQds);
        }
        else{
          f0 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
          runfile_fitted->cd();
          hdQds->Write();
        }
        if(mpv_cosmics > 0 && f0[0]){
          mpv_run_ByLEMs_view0[lem]->Fill(f0[0]/mpv_cosmics);
          if(f0[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
          if(f0[0]/mpv_cosmics > max){max = f0[0]/mpv_cosmics;}
        }
        else if(f0[0] > 0){mpv_run_ByLEMs_view0[lem]->Fill(f0[0]);}
        runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view1").data(), hdQds);
        if(read_fit){
          f1 = ReadFit(hdQds);
        }
        else{
          f1 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
          runfile_fitted->cd();
          hdQds->Write();
        }
        if(mpv_cosmics > 0 && f1[0]){
          mpv_run_ByLEMs_view1[lem]->Fill(f1[0]/mpv_cosmics);
          if(f1[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
          if(f1[0]/mpv_cosmics > max){max = f1[0]/mpv_cosmics;}
        }
        else if(f1[0] > 0){mpv_run_ByLEMs_view1[lem]->Fill(f1[0]);}
        if(mpv_cosmics > 0 && f0[0] > 0 && f1[0] > 0){
          mpv_run_ByLEMs[lem]->Fill((f0[0]+f1[0])/mpv_cosmics);
          if((f0[0]+f1[0])/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        }
        else if(f0[0] > 0 && f1[0] > 0){mpv_run_ByLEMs[lem]->Fill(f0[0]+f1[0]);}
        
        
        runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view0_Dx_Corrected").data(), hdQds);
        if(read_fit){
          f0 = ReadFit(hdQds);
        }
        else{
          f0 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
          runfile_fitted->cd();
          hdQds->Write();
        }
        if(mpv_cosmics > 0 && f0[0]){
          mpv_run_ByLEMs_Dx_Corrected_view0[lem]->Fill(f0[0]/mpv_cosmics);
          if(f0[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
          if(f0[0]/mpv_cosmics > max){max = f0[0]/mpv_cosmics;}
        }
        else if(f0[0] > 0){mpv_run_ByLEMs_Dx_Corrected_view0[lem]->Fill(f0[0]);}
        runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view1_Dx_Corrected").data(), hdQds);
        if(read_fit){
          f1 = ReadFit(hdQds);
        }
        else{
          f1 = fit_dQds(hdQds, false, min_number_of_hits, 0.05, .5);
          runfile_fitted->cd();
          hdQds->Write();
        }
        if(mpv_cosmics > 0 && f1[0]){
          mpv_run_ByLEMs_Dx_Corrected_view1[lem]->Fill(f1[0]/mpv_cosmics);
          if(f1[0]/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
          if(f1[0]/mpv_cosmics > max){max = f1[0]/mpv_cosmics;}
        }
        else if(f1[0] > 0){mpv_run_ByLEMs_Dx_Corrected_view1[lem]->Fill(f1[0]);}
        if(mpv_cosmics > 0 && f0[0] > 0 && f1[0] > 0){
          mpv_run_ByLEMs_Dx_Corrected[lem]->Fill((f0[0]+f1[0])/mpv_cosmics);
          if((f0[0]+f1[0])/mpv_cosmics > gain_max){cout << "WARNING: exceeded max gain of " << gain_max << ". Increase this limit." << endl;}
        }
        else if(f0[0] > 0 && f1[0] > 0){mpv_run_ByLEMs_Dx_Corrected[lem]->Fill(f0[0]+f1[0]);}
      }//for lem
      if(runfile_fitted->IsOpen()){
        runfile_fitted->Close();
      }
      delete runfile_fitted;
      #if verbose
      cout << "    done reading histograms and filling graphs" << endl;
      #endif
      delete hdQds;
      runfile.Close();
      
      #if verbose
      cout << "    next run" << endl;
      #endif
    }//for runs

    if(files_not_found == run_list.size()){
      #if verbose
      cout << "  Not enough files for " << scan_num << "." << endl;
      #endif
      #if verbose
      cout << "    next scan num" << endl;
      #endif
      continue;
    }
    
  
    string outpath = gain_stability_Output;
    check_and_mkdir(outpath);
    outpath = outpath + cut_type;
    check_and_mkdir(outpath);
    outpath = outpath + "/" + to_string(scan_num) + "/";
    check_and_mkdir(outpath);
    string outfile = outpath + "gain.root";
    TFile ofile(outfile.data(), "RECREATE");
    
    #if verbose
    cout << "    Writing file: " << outfile << endl;
    #endif
    
    if( mpv_run->GetEntries() < 2 ){
      #if verbose
      cout << "    No point in MPV/Gain graph " << mpv_run->GetName() << endl;
      #endif
      delete mpv_run; delete mpv_run_view0; delete mpv_run_view1;
      delete mpv_run_Dx_Corrected; delete mpv_run_Dx_Corrected_view0; delete mpv_run_Dx_Corrected_view1;
    }
    else{
      //write the graphs on the root file
      #if verbose
      cout << "    Writing MPV/Gain graph " << mpv_run->GetName() << endl;
      #endif
      
      ofile.cd();
      mpv_run->GetXaxis()->SetTitle("gain");
      mpv_run->GetYaxis()->SetTitle("# runs");
      mpv_run->SetMarkerStyle(29);
      mpv_run->SetMarkerSize(2);
      mpv_run->SetMarkerColor(kRed);
      mpv_run->Write();
      mpv_run->Draw("P");
      gPad->Write(string(string(mpv_run->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run->GetName())+".png").data());
      delete mpv_run;
      delete gPad;
      
      ofile.cd();
      mpv_run_view0->GetXaxis()->SetTitle("gain");
      mpv_run_view0->GetYaxis()->SetTitle("# runs");
      mpv_run_view0->SetMarkerStyle(29);
      mpv_run_view0->SetMarkerSize(2);
      mpv_run_view0->SetMarkerColor(kRed);
      mpv_run_view0->Write();
      mpv_run_view0->Draw("P");
      gPad->Write(string(string(mpv_run_view0->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_view0->GetName())+".png").data());
      delete mpv_run_view0;
      delete gPad;
      
      ofile.cd();
      mpv_run_view1->GetXaxis()->SetTitle("gain");
      mpv_run_view1->GetYaxis()->SetTitle("# runs");
      mpv_run_view1->SetMarkerStyle(29);
      mpv_run_view1->SetMarkerSize(2);
      mpv_run_view1->SetMarkerColor(kRed);
      mpv_run_view1->Write();
      mpv_run_view1->Draw("P");
      gPad->Write(string(string(mpv_run_view1->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_view1->GetName())+".png").data());
      delete mpv_run_view1;
      delete gPad;
      
      ofile.cd();
      mpv_run_Dx_Corrected->GetXaxis()->SetTitle("gain");
      mpv_run_Dx_Corrected->GetYaxis()->SetTitle("# runs");
      mpv_run_Dx_Corrected->SetMarkerStyle(29);
      mpv_run_Dx_Corrected->SetMarkerSize(2);
      mpv_run_Dx_Corrected->SetMarkerColor(kRed);
      mpv_run_Dx_Corrected->Write();
      mpv_run_Dx_Corrected->Draw("P");
      gPad->Write(string(string(mpv_run_Dx_Corrected->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_Dx_Corrected->GetName())+".png").data());
      delete mpv_run_Dx_Corrected;
      delete gPad;
      
      ofile.cd();
      mpv_run_Dx_Corrected_view0->GetXaxis()->SetTitle("gain");
      mpv_run_Dx_Corrected_view0->GetYaxis()->SetTitle("# runs");
      mpv_run_Dx_Corrected_view0->SetMarkerStyle(29);
      mpv_run_Dx_Corrected_view0->SetMarkerSize(2);
      mpv_run_Dx_Corrected_view0->SetMarkerColor(kRed);
      mpv_run_Dx_Corrected_view0->Write();
      mpv_run_Dx_Corrected_view0->Draw("P");
      gPad->Write(string(string(mpv_run_Dx_Corrected_view0->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_Dx_Corrected_view0->GetName())+".png").data());
      delete mpv_run_Dx_Corrected_view0;
      delete gPad;
      
      ofile.cd();
      mpv_run_Dx_Corrected_view1->GetXaxis()->SetTitle("gain");
      mpv_run_Dx_Corrected_view1->GetYaxis()->SetTitle("# runs");
      mpv_run_Dx_Corrected_view1->SetMarkerStyle(29);
      mpv_run_Dx_Corrected_view1->SetMarkerSize(2);
      mpv_run_Dx_Corrected_view1->SetMarkerColor(kRed);
      mpv_run_Dx_Corrected_view1->Write();
      mpv_run_Dx_Corrected_view1->Draw("P");
      gPad->Write(string(string(mpv_run_Dx_Corrected_view1->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_Dx_Corrected_view1->GetName())+".png").data());
      delete mpv_run_Dx_Corrected_view1;
      delete gPad;
    }
    for( auto lem : lems ){
      if( mpv_run_ByLEMs[lem]->GetEntries() < 2 ){
        #if verbose
        cout << "    No point in MPV/Gain graph for LEM " << lem << endl;
        #endif
        delete mpv_run_ByLEMs[lem];
        delete mpv_run_ByLEMs_Dx_Corrected[lem];
        continue;
      }
      #if verbose
      cout << "    Writing MPV/Gain graph for LEM " << lem << endl;
      #endif
      ofile.cd();
      mpv_run_ByLEMs[lem]->GetXaxis()->SetTitle("gain");
      mpv_run_ByLEMs[lem]->GetYaxis()->SetTitle("##runs");
      mpv_run_ByLEMs[lem]->SetMarkerStyle(29);
      mpv_run_ByLEMs[lem]->SetMarkerSize(2);
      mpv_run_ByLEMs[lem]->SetMarkerColor(kRed);
      mpv_run_ByLEMs[lem]->Write();
      mpv_run_ByLEMs[lem]->Draw("P");
      gPad->Write(string(string(mpv_run_ByLEMs[lem]->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_ByLEMs[lem]->GetName())+".png").data());
      delete mpv_run_ByLEMs[lem];
      delete gPad;
      
      ofile.cd();
      mpv_run_ByLEMs_view0[lem]->GetXaxis()->SetTitle("gain");
      mpv_run_ByLEMs_view0[lem]->GetYaxis()->SetTitle("##runs");
      mpv_run_ByLEMs_view0[lem]->SetMarkerStyle(29);
      mpv_run_ByLEMs_view0[lem]->SetMarkerSize(2);
      mpv_run_ByLEMs_view0[lem]->SetMarkerColor(kRed);
      mpv_run_ByLEMs_view0[lem]->Write();
      mpv_run_ByLEMs_view0[lem]->Draw("P");
      gPad->Write(string(string(mpv_run_ByLEMs_view0[lem]->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_ByLEMs_view0[lem]->GetName())+".png").data());
      delete mpv_run_ByLEMs_view0[lem];
      delete gPad;
      
      ofile.cd();
      mpv_run_ByLEMs_view1[lem]->GetXaxis()->SetTitle("gain");
      mpv_run_ByLEMs_view1[lem]->GetYaxis()->SetTitle("##runs");
      mpv_run_ByLEMs_view1[lem]->SetMarkerStyle(29);
      mpv_run_ByLEMs_view1[lem]->SetMarkerSize(2);
      mpv_run_ByLEMs_view1[lem]->SetMarkerColor(kRed);
      mpv_run_ByLEMs_view1[lem]->Write();
      mpv_run_ByLEMs_view1[lem]->Draw("P");
      gPad->Write(string(string(mpv_run_ByLEMs_view1[lem]->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_ByLEMs_view1[lem]->GetName())+".png").data());
      delete mpv_run_ByLEMs_view1[lem];
      delete gPad;
      
      ofile.cd();
      mpv_run_ByLEMs_Dx_Corrected[lem]->GetXaxis()->SetTitle("gain");
      mpv_run_ByLEMs_Dx_Corrected[lem]->GetYaxis()->SetTitle("##runs");
      mpv_run_ByLEMs_Dx_Corrected[lem]->SetMarkerStyle(29);
      mpv_run_ByLEMs_Dx_Corrected[lem]->SetMarkerSize(2);
      mpv_run_ByLEMs_Dx_Corrected[lem]->SetMarkerColor(kRed);
      mpv_run_ByLEMs_Dx_Corrected[lem]->Write();
      mpv_run_ByLEMs_Dx_Corrected[lem]->Draw("P");
      gPad->Write(string(string(mpv_run_ByLEMs_Dx_Corrected[lem]->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_ByLEMs_Dx_Corrected[lem]->GetName())+".png").data());
      delete mpv_run_ByLEMs_Dx_Corrected[lem];
      delete gPad;
      
      ofile.cd();
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->GetXaxis()->SetTitle("gain");
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->GetYaxis()->SetTitle("##runs");
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->SetMarkerStyle(29);
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->SetMarkerSize(2);
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->SetMarkerColor(kRed);
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->Write();
      mpv_run_ByLEMs_Dx_Corrected_view0[lem]->Draw("P");
      gPad->Write(string(string(mpv_run_ByLEMs_Dx_Corrected_view0[lem]->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_ByLEMs_Dx_Corrected_view0[lem]->GetName())+".png").data());
      delete mpv_run_ByLEMs_Dx_Corrected_view0[lem];
      delete gPad;
      
      ofile.cd();
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->GetXaxis()->SetTitle("gain");
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->GetYaxis()->SetTitle("##runs");
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->SetMarkerStyle(29);
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->SetMarkerSize(2);
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->SetMarkerColor(kRed);
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->Write();
      mpv_run_ByLEMs_Dx_Corrected_view1[lem]->Draw("P");
      gPad->Write(string(string(mpv_run_ByLEMs_Dx_Corrected_view1[lem]->GetName())+"_pad").data());
      gPad->SaveAs(string(outpath + string(mpv_run_ByLEMs_Dx_Corrected_view1[lem]->GetName())+".png").data());
      delete mpv_run_ByLEMs_Dx_Corrected_view1[lem];
      delete gPad;
    }//for lem
    ofile.Close();
    
    #if verbose
    cout << "    next scan number" << endl;
    #endif
  }//for scan num
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
