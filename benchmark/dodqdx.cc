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

// vector<int> run_list = {/*1165, 1166, 1167, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1180, 1181, */1182, 1183, 1187, 1188, 1189, 1190, 1191, /*1192, 1193,*/ 1194, 1195, 1196, 1197, /*1198, 1199,*/ 840}
//vector<int> run_list = {1183}
void dodqdx(vector<int> run_list = {/*1165, 1166, 1167, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1180, 1181, */1182, 1183, 1187, 1188, 1189, 1190, 1191, /*1192, 1193,*/ 1194, 1195, 1196, 1197, /*1198, 1199,*/ 840}, string cut_type = "philippe", string v = "July", string m_dQ = "sum", string m_ds = "local", bool save_plots = true){
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
  
  bool recreate_fit_file = true;
  bool update_header = true;
  bool clean = true;
  
  string corrected_or_not = "";
  if(cut_type.find("tg") != string::npos){corrected_or_not = "_Dx_Corrected";}
  string name_to_get = "";
  string cut_type_and_methods = cut_type + "_" + method_ds + "_" + method_dQ;
  
  //****************************************************************************
  
  for(auto run : run_list){
    
    #if verbose
    cout << "  Processing run " << run << endl;
    #endif

    TMyFileHeader header = load_run_header(run,clean);
    if(header.GetRun() == -1){return;}
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
    TF1 func0;
    TF1 func1;
    vector<pair<TH1D,TH1D> > histograms;
    vector<pair<TF1,TF1> > functions;
    vector<pair<double,double> > MPVs;
    vector<pair<double,double> > LMPVs;
    vector<pair<double,double> > LWs;
    vector<pair<double,double> > GWs;
    vector<pair<double,double> > errLMPVs;
    vector<pair<double,double> > errLWs;
    vector<pair<double,double> > errGWs;
    TString FunName = "";
    int nhits_tot = 0;
    int min_number_of_hits = 200;
    vector<double> f = {-1,-1};
    double MPV0, MPV1;
    double LMPV0, LMPV1, err_LMPV0, err_LMPV1;
    double LW0, LW1, err_LW0, err_LW1;
    double GW0, GW1, err_GW0, err_GW1;
    
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
    
    hdQds0.SetTitle("dQds");
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
    
    hdQds1.SetTitle("dQds");
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
      if(func0.GetNpar() > 1){
        MPV0 = func0.GetMaximumX();
        LMPV0 = func0.GetParameter(1);
        LW0 = func0.GetParameter(0);
        GW0 = func0.GetParameter(3);
        err_LMPV0 = func0.GetParError(1);
        err_LW0 = func0.GetParError(0);
        err_GW0 = func0.GetParError(3);
      }
    }
    else{
      func0 = TF1("empty","");
    }
    if( ((TList*)hdQds1.GetListOfFunctions())->GetSize() > 0){
      func1 = *(TF1*)((TList*)hdQds1.GetListOfFunctions())->At(0);
      if(func1.GetNpar() > 1){
        MPV1 = func1.GetMaximumX();
        LMPV1 = func1.GetParameter(1);
        LW1 = func1.GetParameter(0);
        GW1 = func1.GetParameter(3);
        err_LMPV1 = func1.GetParError(1);
        err_LW1 = func1.GetParError(0);
        err_GW1 = func1.GetParError(3);
      }
    } 
    else{
      func1 = TF1("empty","");
    }
      
    histograms.push_back(make_pair(hdQds0,hdQds1) );
    functions.push_back( make_pair(func0,func1) );
    MPVs.push_back( make_pair(MPV0, MPV1) );
    LMPVs.push_back( make_pair(LMPV0, LMPV1) );
    LWs.push_back( make_pair(LW0, LW1) );
    GWs.push_back( make_pair(GW0, GW1) );
    errLMPVs.push_back( make_pair(err_LMPV0, err_LMPV1) );
    errLWs.push_back( make_pair(err_LW0, err_LW1) );
    errGWs.push_back( make_pair(err_GW0, err_GW1) );
    header.SetMPV0(cut_type_and_methods, LMPV0);
    header.SetMPV1(cut_type_and_methods, LMPV1);
    header.SetErrMPV0(cut_type_and_methods, err_LMPV0);
    header.SetErrMPV1(cut_type_and_methods, err_LMPV1);
    header.SetH0(cut_type_and_methods, &hdQds0);
    header.SetH1(cut_type_and_methods, &hdQds1);
    MPV0 = -1; MPV1 = -1;
    LMPV0 = -1; LMPV1 = -1; err_LMPV0 = -1; err_LMPV1 = -1;
    LW0 = -1; LW1 = -1; err_LW0 = -1; err_LW1 = -1;
    GW0 = -1; GW1 = -1; err_GW0 = -1; err_GW1 = -1;
    
    header.ComputeGain(cut_type_and_methods, mpv_cosmics);
    
    #if verbose
    cout << "    Processing " << nhits_tot << " hits for run " << run << endl;
    #endif
    
    for( auto lem : lems ){
    
      name_to_get = "dQds_LEM_"+to_string(lem)+"_view0" + corrected_or_not;
      if(!get_histo_in_inputfile(hdQds0, runfile, name_to_get, read_fit)){return;}
      if(!read_fit && !runfile_fitted->IsOpen()){
        runfile->Close(); delete runfile; runfile = 0;
        runfile = TFile::Open(filename_nonfitted.data(), "READ");
        runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
        if(!get_histo_in_inputfile(hdQds0, runfile, name_to_get, read_fit)){return;}
      }
      
      hdQds0.SetTitle(string("dQds_LEM_"+to_string(lem)).data());
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
      
      hdQds1.SetTitle(string("dQds_LEM_"+to_string(lem)).data());
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
        if(func0.GetNpar() > 1){
          MPV0 = func0.GetMaximumX();
          LMPV0 = func0.GetParameter(1);
          LW0 = func0.GetParameter(0);
          GW0 = func0.GetParameter(3);
          err_LMPV0 = func0.GetParError(1);
          err_LW0 = func0.GetParError(0);
          err_GW0 = func0.GetParError(3);
        }
      }
      else{
        func0 = TF1("empty","");
      }
      if( ((TList*)hdQds1.GetListOfFunctions())->GetSize() > 0){
        func1 = *(TF1*)((TList*)hdQds1.GetListOfFunctions())->At(0);
        if(func1.GetNpar() > 1){
          MPV1 = func1.GetMaximumX();
          LMPV1 = func1.GetParameter(1);
          LW1 = func1.GetParameter(0);
          GW1 = func1.GetParameter(3);
          err_LMPV1 = func1.GetParError(1);
          err_LW1 = func1.GetParError(0);
          err_GW1 = func1.GetParError(3);
        }
      }
      else{
        func1 = TF1("empty","");
      }
      histograms.push_back(make_pair(hdQds0,hdQds1) );
      functions.push_back( make_pair(func0,func1) );
      MPVs.push_back( make_pair(MPV0, MPV1) );
      LMPVs.push_back( make_pair(LMPV0, LMPV1) );
      LWs.push_back( make_pair(LW0, LW1) );
      GWs.push_back( make_pair(GW0, GW1) );
      errLMPVs.push_back( make_pair(err_LMPV0, err_LMPV1) );
      errLWs.push_back( make_pair(err_LW0, err_LW1) );
      errGWs.push_back( make_pair(err_GW0, err_GW1) );
      header.SetMPVLEM(cut_type_and_methods, lem, make_pair(LMPV0,LMPV1));
      header.SetErrMPVLEM(cut_type_and_methods, lem, make_pair(err_LMPV0,err_LMPV1));
      header.SetHLEM(cut_type_and_methods, lem, make_pair(&hdQds0,&hdQds1));
      MPV0 = -1; MPV1 = -1; err_MPV0 = -1; err_MPV1 = -1;
      LMPV0 = -1; LMPV1 = -1; err_LMPV0 = -1; err_LMPV1 = -1;
      LW0 = -1; LW1 = -1; err_LW0 = -1; err_LW1 = -1;
      GW0 = -1; GW1 = -1; err_GW0 = -1; err_GW1 = -1;
      
      header.ComputeGain(cut_type_and_methods, mpv_cosmics,lem);
    
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
    if(update_header){
      save_run_header(header);
    }
    
    if(!save_plots){return;}
    
    
    string outpath = dQds_Output + cut_type + "/plots/";
    check_and_mkdir(outpath);
    string outfile = outpath + to_string(run) + "_dQds";
    for(int i = 0; i < histograms.size(); i++){
      TCanvas mycan("","",1000,620);
      TLegend leg(0.75,0.7,0.95,0.95);
      histograms[i].second.SetLineColor(kBlue);
      functions[i].second.SetLineColor(kBlue);
      string title0 = "#splitline{v0 Langau MPV: " + to_string_with_precision(MPVs[i].first,2) + "}{#splitline{Landau MPV: " + to_string_with_precision(LMPVs[i].first,2) + "+/-" + to_string_with_precision(errLMPVs[i].first,2) + "}{#splitline{LWidth: " + to_string_with_precision(LWs[i].first,2) + "+/-" + to_string_with_precision(errLWs[i].first,2) + "}{GWidth: " + to_string_with_precision(GWs[i].first,2) + "+/-" + to_string_with_precision(errGWs[i].first,2) + "}}}";
      string title1 = "#splitline{v1 Langau MPV: " + to_string_with_precision(MPVs[i].second,2) + "}{#splitline{Landau MPV: " + to_string_with_precision(LMPVs[i].second,2) + "+/-" + to_string_with_precision(errLMPVs[i].second,2) + "}{#splitline{LWidth: " + to_string_with_precision(LWs[i].second,2) + "+/-" + to_string_with_precision(errLWs[i].second,2) + "}{GWidth: " + to_string_with_precision(GWs[i].second,2) + "+/-" + to_string_with_precision(errGWs[i].second,2) + "}}}";
      double max = histograms[i].second.GetMaximum();
      if(histograms[i].first.GetMaximum() > max){max = histograms[i].first.GetMaximum();}
      histograms[i].first.SetLineColor(kRed);
      functions[i].first.SetLineColor(kRed);
      histograms[i].second.SetMaximum(max*1.1);
      histograms[i].second.Draw();
      histograms[i].first.Draw("SAME");
      leg.AddEntry(&histograms[i].first,title0.data());
      leg.AddEntry(&histograms[i].second,title1.data());
      leg.Draw();
      if(string(functions[i].second.GetName()) != "empty"){functions[i].second.Draw("SAME");}
      if(string(functions[i].first.GetName()) != "empty"){functions[i].first.Draw("SAME");}
      if(i == 0){
        mycan.Print(string(outfile+".pdf(").data(),"pdf");
      }
      else if(i == histograms.size()-1){
        mycan.Print(string(outfile+".pdf)").data(),"pdf");
      }
      else{
        mycan.Print(string(outfile+".pdf").data(),"pdf");
      }
    }
    histograms.clear();
    functions.clear();
    MPVs.clear();
    LMPVs.clear();
    LWs.clear();
    GWs.clear();
    errLMPVs.clear();
    errLWs.clear();
    errGWs.clear();
    
    #if verbose
    cout << "    next run" << endl;
    #endif
  }//for runs
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
