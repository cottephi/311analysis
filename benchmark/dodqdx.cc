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


void dodqdx(vector<int> run_list = {}, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D", bool save_plots = true){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  if(!load_run_lists()){return;}
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);
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
  
  bool recreate_fit_file = true;
  
  string corrected_or_not = "";
  if(cut_type.find("tg") != string::npos){corrected_or_not = "_Dx_Corrected";}
  string name_to_get = "";
  
  //****************************************************************************
  load_cosmics();
  
  for(auto run : run_list){

    string filename_nonfitted = dQds_Output+cut_type+"/"+to_string(run)+".root";
    string filename_fitted = dQds_Output+cut_type+"/"+to_string(run)+"_fitted.root";
    TFile *runfile_fitted = new TFile();
    bool read_fit = false;
    string ifilename = "";
    int files_not_found = 0;
    int i_read_fit = read_or_do_fit(filename_fitted, filename_nonfitted, recreate_fit_file, ifilename, files_not_found);
    if(i_read_fit == 3){continue;}
    if(i_read_fit == 1){read_fit = true;}
    TFile *runfile = TFile::Open(ifilename.data(),"READ");
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

    TH1D hdQds0;
    TH1D hdQds1;
    vector<pair<TH1D,TH1D> > histograms;
    vector<pair<TF1,TF1> > functions;
    TString FunName = "";
    int nhits_tot = 0;
    int min_number_of_hits = 200;
    vector<double> f = {-1,-1};
    
    #if verbose 
    cout << "  Reading histograms..." << endl;
    #endif
    
    name_to_get = "dQds_view0" + corrected_or_not;
    if(!get_histo_in_inputfile(hdQds0, runfile, name_to_get, read_fit)){return;}
    if(!read_fit && !runfile_fitted->IsOpen()){
      runfile->Close(); delete runfile; runfile = 0;
      runfile = TFile::Open(filename_nonfitted.data(), "READ");
      runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
      if(!get_histo_in_inputfile(hdQds0, runfile, name_to_get, read_fit)){return;}
    }
    
    nhits_tot += hdQds0.GetEntries();
    if(read_fit){
      f = ReadFit(&hdQds0);
//            if(runs_and_fields[run]["Extraction"] < 1800){}
    }
    else{
      f = fit_dQds(&hdQds0, false, min_number_of_hits, 0.05, 0.5);
      runfile_fitted->cd();
      hdQds0.Write();
    }
    
    name_to_get = "dQds_view1" + corrected_or_not;
    if(!get_histo_in_inputfile(hdQds1, runfile, name_to_get, read_fit)){return;}
    if(!read_fit && !runfile_fitted->IsOpen()){
      runfile->Close(); delete runfile; runfile = 0;
      runfile = TFile::Open(filename_nonfitted.data(), "READ");
      runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
      if(!get_histo_in_inputfile(hdQds1, runfile, name_to_get, read_fit)){return;}
    }
    
    nhits_tot += hdQds1.GetEntries();
    if(read_fit){
      f = ReadFit(&hdQds1);
    }
    else{
      f = fit_dQds(&hdQds1, false, min_number_of_hits, 0.05, 0.5);
      runfile_fitted->cd();
      hdQds1.Write();
    }
    histograms.push_back(make_pair(hdQds0,hdQds1) );
    functions.push_back( make_pair(*(TF1*)((TList*)hdQds0.GetListOfFunctions())->At(0), *(TF1*)((TList*)hdQds1.GetListOfFunctions())->At(0)) );
    
    if( nhits_tot < min_number_of_hits ){
      #if verbose
      cout << "    Insuficient number of hits (" << nhits_tot << ") for run " << run << ". Skipping this run. \n  next run" << endl;
      #endif
      if(runfile_fitted->IsOpen()){
        runfile_fitted->Close();
      }
      delete runfile_fitted;
      runfile->Close();
      delete runfile;
      continue;
    }
    else{
      #if verbose
      cout << "    Processing " << nhits_tot << " hits for run " << run << endl;
      #endif
    }
    
    for( auto lem : lems ){
    
      name_to_get = "dQds_LEM_"+to_string(lem)+"_view0" + corrected_or_not;
      if(!get_histo_in_inputfile(hdQds0, runfile, name_to_get, read_fit)){return;}
      if(!read_fit && !runfile_fitted->IsOpen()){
        runfile->Close(); delete runfile; runfile = 0;
        runfile = TFile::Open(filename_nonfitted.data(), "READ");
        runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
        if(!get_histo_in_inputfile(hdQds0, runfile, name_to_get, read_fit)){return;}
      }
      
      if(read_fit){
        f = ReadFit(&hdQds0);
      }
      else{
        f = fit_dQds(&hdQds0, false, min_number_of_hits, 0.05, 0.5);
        runfile_fitted->cd();
        hdQds0.Write();
      }
    
      name_to_get = "dQds_LEM_"+to_string(lem)+"_view1" + corrected_or_not;
      if(!get_histo_in_inputfile(hdQds1, runfile, name_to_get, read_fit)){return;}
      if(!read_fit && !runfile_fitted->IsOpen()){
        runfile->Close(); delete runfile; runfile = 0;
        runfile = TFile::Open(filename_nonfitted.data(), "READ");
        runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
        if(!get_histo_in_inputfile(hdQds1, runfile, name_to_get, read_fit)){return;}
      }
      
      if(read_fit){
        f = ReadFit(&hdQds1);
      }
      else{
        f = fit_dQds(&hdQds1, false, min_number_of_hits, 0.05, 0.5);
        runfile_fitted->cd();
        hdQds1.Write();
      }
      histograms.push_back(make_pair(hdQds0,hdQds1) );
      functions.push_back( make_pair(*(TF1*)((TList*)hdQds0.GetListOfFunctions())->At(0), *(TF1*)((TList*)hdQds1.GetListOfFunctions())->At(0)) );
    }//for lem
    if(runfile_fitted->IsOpen()){
      runfile_fitted->Close();
    }
    delete runfile_fitted;
    #if verbose
    cout << "    done reading histograms and filling graphs" << endl;
    #endif
    runfile->Close();
    delete runfile;
    
    if(!save_plots){return;}
    
    
    string outpath = dQds_Output + cut_type + "/plots/" + to_string(run) + "/";
    check_and_mkdir(outpath);
    for(int i = 0; i < histograms.size(); i++){
      string outfile = outpath + string(histograms[i].first.GetName()).erase(string(histograms[i].first.GetName()).find("_view"),6);
      histograms[i].second.SetLineColor(kRed);
      functions[i].second.SetLineColor(kRed);
      histograms[i].second.SetTitle("view0 blue, view1 red;fC/cm");
      double max = histograms[i].second.GetMaximum();
      if(histograms[i].first.GetMaximum() > max){max = histograms[i].first.GetMaximum();}
      functions[i].first.SetLineColor(kBlue);
      histograms[i].second.SetMaximum(max*1.1);
      histograms[i].second.Draw();
      functions[i].second.Draw("SAME");
      histograms[i].first.Draw("SAME");
      functions[i].first.Draw("SAME");
      gPad->SaveAs(string(outfile+".png").data());
      delete gPad;
    }
    
    #if verbose
    cout << "    next run" << endl;
    #endif
  }//for runs
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
