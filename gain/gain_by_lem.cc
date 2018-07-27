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

pair<int,int> get_lem_bin(int lem){
  if(lem == 2){return make_pair(0.5,1.5);}
  if(lem == 4){return make_pair(0.5,0.5);}
  if(lem == 5){return make_pair(1.5,1.5);}
  if(lem == 6){return make_pair(2.5,1.5);}
  if(lem == 7){return make_pair(1.5,0.5);}
  if(lem == 8){return make_pair(2.5,0.5);}
  if(lem == 9){return make_pair(3.5,1.5);}
  if(lem == 11){return make_pair(3.5,0.5);}
}

vector<TPaveText*> draw_lem_legend(){
  vector<TPaveText*> lem_legend = {};
  lem_legend.push_back(new TPaveText(0.45,1.45, 0.55,1.55));
  lem_legend.push_back(new TPaveText(0.45,0.45, 0.55,0.55));
  lem_legend.push_back(new TPaveText(1.45,1.45, 1.55,1.55));
  lem_legend.push_back(new TPaveText(2.45,1.45, 2.55,1.55));
  lem_legend.push_back(new TPaveText(1.45,0.45, 1.55,0.55));
  lem_legend.push_back(new TPaveText(2.45,0.45, 2.55,0.55));
  lem_legend.push_back(new TPaveText(3.45,1.45, 3.55,1.55));
  lem_legend.push_back(new TPaveText(3.45,0.45, 3.55,0.55));
  for( int lem = 0; lem < lem_legend.size(); lem++ ){
    lem_legend[lem]->AddText(to_string(lem).data());
  }
  return lem_legend;
}

void gain_by_lem(int run = 840, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  gErrorIgnoreLevel = kError;
//  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);
  
  vector<TPaveText*> lem_legend = draw_lem_legend();
  
//  gStyle->SetPalette(kColorPrintableOnGrey);

  load_cosmics();
  
  TH2D mpv_ByLEMs("mpv_ByLEMs","mpv_ByLEMs",4,0,4,2,0,2);
  TH2D mpv_ByLEMs_view0("mpv_ByLEMs_view0","mpv_ByLEMs_view0",4,0,4,2,0,2);
  TH2D mpv_ByLEMs_view1("mpv_ByLEMs_view1","mpv_ByLEMs_view1",4,0,4,2,0,2);
  
  TH1D distri_ByLEMs("distri_ByLEMs","distri_ByLEMs",100,0,5);
  TH1D distri_ByLEMs_view0("distri_ByLEMs_view0","distri_ByLEMs_view0",100,0,2.5);
  TH1D distri_ByLEMs_view1("distri_ByLEMs_view1","distri_ByLEMs_view1",100,0,2.5);
  
  //Minimum number of hits to have correct statistics
  int min_number_of_hits = 200;


  bool recreate_fit_file = false;
  string filename_nonfitted = dQds_Output+cut_type+"/"+to_string(run)+".root";
  string filename_fitted = dQds_Output+cut_type+"/"+to_string(run)+"_fitted.root";
  TFile *runfile_fitted = new TFile();
  bool read_fit = false;
  string ifilename = "";
  int files_not_found = 0;
  int i_read_fit = read_or_do_fit(filename_fitted, filename_nonfitted, recreate_fit_file, ifilename, files_not_found);
  if(i_read_fit == 3){return;}
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
  
  TH1D hdQds;
  TString FunName = "";
  int nhits_tot = 0;
  vector<double> f0 = {-1,-1};
  vector<double> f1 = {-1,-1};
  string name_to_get = "";
  string corrected_or_not = "";
  if(cut_type.find("throughgoing") != string::npos){corrected_or_not = "_Dx_Corrected";}

  #if verbose 
  cout << "  Reading histograms..." << endl;
  #endif
  for( auto lem : lems){
    pair<int,int> bins = get_lem_bin(lem);
    
    name_to_get = string("dQds_LEM_"+to_string(lem)+"_view0" + corrected_or_not);
    if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
    if(!read_fit && !runfile_fitted->IsOpen()){
      runfile->Close(); delete runfile; runfile = 0;
      runfile = TFile::Open(filename_nonfitted.data(), "READ");
      runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
      if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
    }
    
    if(read_fit){
      f0 = ReadFit(&hdQds);
    }
    else{
      f0 = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, .5);
      runfile_fitted->cd();
      hdQds.Write();
    }
    if(mpv_cosmics > 0 && f0[0] > 0){
      mpv_ByLEMs_view0.Fill(bins.first,bins.second,f0[0]/mpv_cosmics);
      distri_ByLEMs_view0.Fill(f0[0]/mpv_cosmics);
    }
    else if(f0[0] > 0){
      mpv_ByLEMs_view0.Fill(bins.first,bins.second,f0[0]);
      distri_ByLEMs_view0.Fill(f0[0]);
    }
    
    name_to_get = string("dQds_LEM_"+to_string(lem)+"_view1" + corrected_or_not);
    if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
    if(!read_fit && !runfile_fitted->IsOpen()){
      runfile->Close(); delete runfile; runfile = 0;
      runfile = TFile::Open(filename_nonfitted.data(), "READ");
      runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
      if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
    }
    
    if(read_fit){
      f1 = ReadFit(&hdQds);
    }
    else{
      f1 = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, .5);
      runfile_fitted->cd();
      hdQds.Write();
    }
    if(mpv_cosmics > 0 && f1[0] > 0){
      mpv_ByLEMs_view1.Fill(bins.first,bins.second,f1[0]/mpv_cosmics);
      distri_ByLEMs_view1.Fill(f1[0]/mpv_cosmics);
    }
    else if(f1[0] > 0){
      mpv_ByLEMs_view1.Fill(bins.first,bins.second,f1[0]);
      distri_ByLEMs_view1.Fill(f1[0]);
    }
    
    if(mpv_cosmics > 0 && f0[0] > 0 && f1[0] > 0){
      mpv_ByLEMs.Fill(bins.first,bins.second,(f0[0]+f1[0])/mpv_cosmics);
      distri_ByLEMs.Fill((f0[0]+f1[0])/mpv_cosmics);
    }
    else if(f0[0] > 0 && f1[0] > 0){
      mpv_ByLEMs.Fill(bins.first,bins.second,f0[0]+f1[0]);
      distri_ByLEMs.Fill(f0[0]+f1[0]);
    }
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
        
  
  string outpath = gain_by_lem_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type;
  check_and_mkdir(outpath);
  outpath = outpath+"/"+to_string(run) + "/";
  check_and_mkdir(outpath);
  string outfile = outpath + "gain.root";
  TFile ofile(outfile.data(), "RECREATE");
  #if verbose
  cout << "    Writing file: " << string(ofile.GetName()).data() << endl;
  #endif
  #if verbose
  cout << "    Writing MPV/Gain graph " << mpv_ByLEMs.GetName() << endl;
  #endif
  
  gStyle->SetOptStat(0);
  ofile.cd();
  mpv_ByLEMs.SetMinimum(0);
  mpv_ByLEMs.SetMaximum(5);
  mpv_ByLEMs.Draw("ACOLZ");
  for(auto leg : lem_legend){leg->Draw();}
  mpv_ByLEMs.Write();
  gPad->SetName(string(string(mpv_ByLEMs.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(mpv_ByLEMs.GetName()) + ".png").data());

  gStyle->SetOptStat(1111);
  ofile.cd();
  distri_ByLEMs.Draw();
  distri_ByLEMs.Write();
  gPad->SetName(string(string(distri_ByLEMs.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(distri_ByLEMs.GetName()) + ".png").data());
  
  gStyle->SetOptStat(0);
  ofile.cd();
  mpv_ByLEMs_view0.SetMinimum(0);
  mpv_ByLEMs_view0.SetMaximum(5);
  mpv_ByLEMs_view0.Draw("ACOLZ");
  for(auto leg : lem_legend){leg->Draw();}
  mpv_ByLEMs_view0.Write();
  gPad->SetName(string(string(mpv_ByLEMs_view0.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(mpv_ByLEMs_view0.GetName()) + ".png").data());

  gStyle->SetOptStat(1111);
  ofile.cd();
  distri_ByLEMs_view0.Draw();
  distri_ByLEMs_view0.Write();
  gPad->SetName(string(string(distri_ByLEMs_view0.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(distri_ByLEMs_view0.GetName()) + ".png").data());
  
  
  gStyle->SetOptStat(0);
  ofile.cd();
  mpv_ByLEMs_view1.SetMinimum(0);
  mpv_ByLEMs_view1.SetMaximum(5);
  mpv_ByLEMs_view1.Draw("ACOLZ");
  for(auto leg : lem_legend){leg->Draw();}
  mpv_ByLEMs_view1.Write();
  gPad->SetName(string(string(mpv_ByLEMs_view1.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(mpv_ByLEMs_view1.GetName()) + ".png").data());

  gStyle->SetOptStat(1111);
  ofile.cd();
  distri_ByLEMs_view1.Draw();
  distri_ByLEMs_view1.Write();
  gPad->SetName(string(string(distri_ByLEMs_view1.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(distri_ByLEMs_view1.GetName()) + ".png").data());
  
  delete gPad;
  ofile.Close();
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
