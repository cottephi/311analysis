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


void YZ_dodqdx(vector<int> run_list = {}, string cut_type = "Ds", string v = "July", string m_dQ = "sum", string m_ds = "local", bool save_plots = true){
  TH1::AddDirectory(kFALSE);
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(v)){return;}
  if(!load_run_lists()){return;}
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  
//  gStyle->SetOptFit(0);
  if (run_list.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path_wa105_311data << "*..." << endl;
    #endif
    string wildcard_path = path_wa105_311data + "*";
    for( auto irun : glob(wildcard_path) ){run_list.push_back(atoi(irun.data()));}
  }
  if(!load_cosmics()){return;}
  
  bool recreate_fit_file = false;
  bool update_header = true;
  
  int nbinsY = (int)(2*lem_size/dy);
  int nbinsZ = (int)(6*lem_size/dz);
  vector<int> begins = {-(int)(lem_size/dy), 0};
  vector<int> ends = {(int)(lem_size/dy), (int)(6*lem_size/dz)};
  
  string corrected_or_not = "";
  if(cut_type.find("tg") != string::npos){corrected_or_not = "_Dx_Corrected";}
  string name_to_get = "";
  string cut_type_and_methods = cut_type + "_" + method_ds + "_" + method_dQ;
  
  //****************************************************************************
  
  for(auto run : run_list){

    TMyFileHeader header = load_run_header(run);
    if(header.GetRun() == -1){return;}
    string filename_nonfitted = dQds_YZ_Output+cut_type+"/"+to_string(run)+".root";
    string filename_fitted = dQds_YZ_Output+cut_type+"/"+to_string(run)+"_fitted.root";
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
    TF1 func0;
    TF1 func1;
    vector<pair<TH1D,TH1D> > histograms;
    vector<pair<TF1,TF1> > functions;
    vector<pair<double,double> > MPVs;
    TString FunName = "";
    int min_number_of_hits = 200;
    vector<double> f = {-1,-1};
    double MPV0, MPV1;
    
    #if verbose 
    cout << "  Reading histograms..." << endl;
    #endif
    
    for( int Y = begins[0]; Y < ends[0]; Y++ ){
      double y = (Y+0.5)*dy;
      for( int Z = begins[1]; Z < ends[1]; Z++ ){
        int nhits_tot = 0;
        double z = (Z+0.5)*dz;
        string sign = "";
        if(Y < 0){sign = "m";}
        name_to_get = string("dQds_view0_" + corrected_or_not + "dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z));
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
        
        name_to_get = string("dQds_view1_" + corrected_or_not + "dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z));
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
        
        if( ((TList*)hdQds0.GetListOfFunctions())->GetSize() > 0){
          func0 = *(TF1*)((TList*)hdQds0.GetListOfFunctions())->At(0);
          MPV0 = func0.GetMaximumX();
        }
        else{
          func0 = TF1("empty","");
          MPV0 = -1;
        }
        if( ((TList*)hdQds1.GetListOfFunctions())->GetSize() > 0){
          func1 = *(TF1*)((TList*)hdQds1.GetListOfFunctions())->At(0);
          MPV1 = func1.GetMaximumX();
        } 
        else{
          func1 = TF1("empty","");
          MPV1 = -1;
        }
        if(nhits_tot < min_number_of_hits){continue;}
        histograms.push_back(make_pair(hdQds0,hdQds1) );
        functions.push_back( make_pair(func0,func1) );
        MPVs.push_back( make_pair(MPV0, MPV1) );
        header.SetYZMPV(cut_type_and_methods, make_pair(Y,Z), make_pair(MPV0,MPV1));
        header.SetYZH(cut_type_and_methods, make_pair(Y,Z), make_pair(&hdQds0,&hdQds1));
        header.ComputeYZGain(cut_type_and_methods, mpv_cosmics, Y, Z);
        
      }//for Z
    }//for Y
    if(runfile_fitted->IsOpen()){
      runfile_fitted->Close();
    }
    delete runfile_fitted;
    #if verbose
    cout << "    done reading histograms and filling graphs" << endl;
    #endif
    runfile->Close();
    delete runfile;
    if(update_header){
      save_run_header(header);
    }
    if(!save_plots){return;}


    string outpath = dQds_YZ_Output + cut_type + "/plots/" + to_string(run) + "/";
    check_and_mkdir(outpath);
    for(int i = 0; i < histograms.size(); i++){
      TCanvas mycan("","",1000,620);
      string outfile = outpath + string(histograms[i].first.GetName()).erase(string(histograms[i].first.GetName()).find("_view"),6);
      histograms[i].second.SetLineColor(kRed);
      functions[i].second.SetLineColor(kRed);
      histograms[i].second.SetTitle(string("view0 blue MPV "+to_string(MPVs[i].first)+" " + to_string((int)histograms[i].first.GetEntries()) + " hits, view1 "+to_string(MPVs[i].second)+" " + to_string((int)histograms[i].second.GetEntries()) + " hits red;fC/cm").data());
      double max = histograms[i].second.GetMaximum();
      if(histograms[i].first.GetMaximum() > max){max = histograms[i].first.GetMaximum();}
      functions[i].first.SetLineColor(kBlue);
      histograms[i].second.SetMaximum(max*1.1);
      histograms[i].second.Draw();
      if(string(functions[i].second.GetName()) != "empty"){functions[i].second.Draw("SAME");}
      if(string(functions[i].first.GetName()) != "empty"){functions[i].first.Draw("SAME");}
      histograms[i].first.Draw("SAME");
      mycan.SaveAs(string(outfile+".png").data());
    }
    histograms.clear();
    functions.clear();
    MPVs.clear();
    
    #if verbose
    cout << "    next run" << endl;
    #endif
  }//for runs
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
