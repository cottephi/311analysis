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


void gain_ana_2D(string cut_type = "Ds", string version = "Feb", string m_dQ = "sum", string m_ds = "local"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  if(!load_run_lists()){return;}
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  
  gStyle->SetPalette(kColorPrintableOnGrey);
  bool recreate_fit_file = false;
  
  string corrected_or_not = "";
  if(cut_type.find("tg") != string::npos){corrected_or_not = "_Dx_Corrected";}
  string name_to_get = "";
  
  //****************************************************************************
  //variables definition here
  map<string, vector<int> > scans;
  
//Generated with python macro (all runs), 50% drift tolerence
//  scans["Amplification_Extraction"] = {2};
  scans["Amplification_Extraction"] = {1,2,3,4,5,6};
  scans["Amplification_Induction"] = {1,2,3,4};
//  scans["Induction_Extraction"] = {1,2,3,4,5,6};
  
  //system(string("rm "+GainAna_Output+"*.root").data());
  load_cosmics();
  
  string outpath = GainAna_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type;
  check_and_mkdir(outpath);
  
  for( auto scan_it : scans ){
  
    string scan_type = scan_it.first;
    outpath = GainAna_Output + cut_type + "/" + scan_type;
    check_and_mkdir(outpath);
    double x0,x1,y0,y1;
    if(scan_type == "Amplification_Extraction"){
      x0 = 21; x1 = 31; y0 = 1; y1 = 3;
    }
    else if(scan_type == "Amplification_Induction"){
      x0 = 21; x1 = 31; y0 = 1; y1 = 5;
    }
    else if(scan_type == "Induction_Extraction"){
      x0 = 1; x1 = 5; y0 = 1; y1 = 3;
    }
    
    string field_type = scan_type;
    string unit = "kV/cm";
      
    int empty_scan_num = 0;
    for( auto scan_num : scan_it.second ){

      //multigraph for mpv, containing view 0 and view 1
      TGraph2D* mpv_field = new TGraph2D();
      TGraph2D* mpv_field_view0 = new TGraph2D();
      TGraph2D* mpv_field_view1 = new TGraph2D();
      TH2D* h_mpv_field = new TH2D(string("h_gain_"+scan_type+"_"+to_string(scan_num)).data(),to_string(scan_num).data(),10,x0,x1,10,y0,y1);
      TH2D* h_mpv_field_view0 = new TH2D(string("h_gain_"+scan_type+"_"+to_string(scan_num)+"_view0").data(),to_string(scan_num).data(),10,x0,x1,10,y0,y1);
      TH2D* h_mpv_field_view1 = new TH2D(string("h_gain_"+scan_type+"_"+to_string(scan_num)+"_view1").data(),to_string(scan_num).data(),10,x0,x1,10,y0,y1);
      mpv_field->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)).data());
      mpv_field->SetTitle(to_string(scan_num).data());
      mpv_field_view0->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_view0").data());
      mpv_field_view0->SetTitle(to_string(scan_num).data());
      mpv_field_view1->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_view1").data());
      mpv_field_view1->SetTitle(to_string(scan_num).data());
      map<int, TGraph2D*> mpv_field_ByLEMs;
      map<int, TGraph2D*> mpv_field_ByLEMs_view0;
      map<int, TGraph2D*> mpv_field_ByLEMs_view1;
      map<int, TH2D*> h_mpv_field_ByLEMs;
      map<int, TH2D*> h_mpv_field_ByLEMs_view0;
      map<int, TH2D*> h_mpv_field_ByLEMs_view1;
      for( auto lem : lems ){
        mpv_field_ByLEMs[lem] = new TGraph2D();
        mpv_field_ByLEMs_view0[lem] = new TGraph2D();
        mpv_field_ByLEMs_view1[lem] = new TGraph2D();
        h_mpv_field_ByLEMs[lem] = new TH2D(string("h_gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)).data(),to_string(scan_num).data(),10,x0,x1,10,y0,y1);
        h_mpv_field_ByLEMs_view0[lem] = new TH2D(string("h_gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_view0").data(),to_string(scan_num).data(),10,x0,x1,10,y0,y1);
        h_mpv_field_ByLEMs_view1[lem] = new TH2D(string("h_gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_view1").data(),to_string(scan_num).data(),10,x0,x1,10,y0,y1);
        mpv_field_ByLEMs[lem]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)).data());
        mpv_field_ByLEMs_view0[lem]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_view0").data());
        mpv_field_ByLEMs_view1[lem]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_view1").data());
        mpv_field_ByLEMs[lem]->SetTitle(to_string(scan_num).data());
        mpv_field_ByLEMs_view0[lem]->SetTitle(to_string(scan_num).data());
        mpv_field_ByLEMs_view1[lem]->SetTitle(to_string(scan_num).data());
      }
      #if verbose 
      cout << "  Initialising graphs..." << endl;
      #endif
      TFile *ofile = new TFile();
      outpath = GainAna_Output+cut_type+"/"+scan_type+"/"+to_string(scan_num)+"/";
      check_and_mkdir(outpath);
      outpath = outpath+"gain.root";
      ofile = TFile::Open(outpath.data(), "RECREATE");

      //Minimum number of hits to have correct statistics
      int min_number_of_hits = 200;

      vector<int> run_list;
      vector<pair<float,float> > field;
      vector<string> const_fields_names;
      vector<float> const_fields_values;

      //12 drift field variation tolerance
      // Drift :
      
      if(!load_runs_2D(scan_type, scan_num, run_list, field, const_fields_names, const_fields_values)){continue;}
      
      #if verbose
      cout << "  Doing " << scan_type << " scan number " << scan_num << endl;
      #endif
      
      //****************************************************************************
      //define an output file

      const int n_runs = run_list.size();
      
      //****************************************************************************
      //create efiled2run map

      map<pair<float, float>, vector<int>> field2run;

      for( int ii=0; ii < n_runs; ii++ ){
        field2run[ field[ii] ].push_back( run_list[ii] );
      }

      //****************************************************************************
      //read files and fit the distributions
      
      int empty_fields = 0;
      
      for( auto field_it : field2run ){
      
        int found_run_for_field = 0;

        #if verbose
        cout << "    Processing fields: " << field_it.first.first << " and " << field_it.first.second << endl;
        #endif

        int files_not_found = 0;
        for(auto run : field_it.second){
          cout << "      procesing run " << run << endl;
          if(!load_rho_run(run)){continue;}
          
          #if verbose
          cout << "      Density correction: " << gain_correction_for_rho(rho_run,runs_and_fields[run]["Amplification"]) << endl;
          #endif

          string filename_nonfitted = dQds_Output+cut_type+"/"+to_string(run)+".root";
          string filename_fitted = dQds_Output+cut_type+"/"+to_string(run)+"_fitted.root";
          TFile *runfile_fitted = new TFile();
          bool read_fit = false;
          string ifilename = "";
          int i_read_fit = read_or_do_fit(filename_fitted, filename_nonfitted, recreate_fit_file, ifilename, files_not_found);
          if(i_read_fit == 3){continue;}
          if(i_read_fit == 1){read_fit = true;}
          TFile *runfile = TFile::Open(ifilename.data(),"READ");
          if( !read_fit ){
            #if verbose
            cout << "      Recording fitted histograms in " << filename_fitted << endl;
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
          double f0 = -1;
          double f1 = -1;
          
          #if verbose 
          cout << "      Reading histograms..." << endl;
          #endif
          
          name_to_get = "dQds_view0" + corrected_or_not;
          if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
          if(!read_fit && !runfile_fitted->IsOpen()){
            runfile->Close(); delete runfile; runfile = 0;
            runfile = TFile::Open(filename_nonfitted.data(), "READ");
            runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
            if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
          }
          
          nhits_tot += hdQds.GetEntries();
          if(read_fit){
            f0 = ReadFit(&hdQds)[0];
          }
          else{
            f0 = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, .5)[0];
            runfile_fitted->cd();
            hdQds.Write();
          }
          if(f0 < 0 and force_mpv.find(run) != force_mpv.end()){
            #if verbose
            cout << "Forcing MPV value..." << endl;
            #endif
            f0 = force_mpv[run].first;
          }
          if(f0 > 0){
            fill_2d_graph(mpv_field_view0, field_it.first.first, field_it.first.second, f0);
            fill_2d_hist(h_mpv_field_view0, field_it.first.first, field_it.first.second, f0);
          }
          
          name_to_get = "dQds_view1" + corrected_or_not;
          if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
          if(!read_fit && !runfile_fitted->IsOpen()){
            runfile->Close(); delete runfile; runfile = 0;
            runfile = TFile::Open(filename_nonfitted.data(), "READ");
            runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
            if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
          }
          
          nhits_tot += hdQds.GetEntries();
          if(read_fit){
            f1 = ReadFit(&hdQds)[0];
          }
          else{
            f1 = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, .5)[0];
            runfile_fitted->cd();
            hdQds.Write();
          }
          if(f1 < 0 and force_mpv.find(run) != force_mpv.end()){
            #if verbose
            cout << "Forcing MPV value..." << endl;
            #endif
            f1 = force_mpv[run].second;
          }
          if(f1 > 0){
            fill_2d_graph(mpv_field_view1, field_it.first.first, field_it.first.second, f1);
            fill_2d_hist(h_mpv_field_view1, field_it.first.first, field_it.first.second, f1);
          }
          if(f0 > 0 && f1 > 0){
            fill_2d_graph(mpv_field, field_it.first.first, field_it.first.second, f0+f1);
            fill_2d_hist(h_mpv_field, field_it.first.first, field_it.first.second, f0+f1);
          }
          if( nhits_tot < min_number_of_hits ){
            #if verbose
            cout << "      Insuficient number of hits (" << nhits_tot << ") for run " << run << ". Skipping this run. \nnext run" << endl;
            #endif
            runfile->Close();
            delete runfile;
            continue;
          }
          else{
            #if verbose
            cout << "      Processing " << nhits_tot << " hits for run " << run << endl;
            #endif
          }
          
          
          for( auto lem : lems ){
            name_to_get = "dQds_LEM_"+to_string(lem)+"_view0" + corrected_or_not;
            if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
            if(!read_fit && !runfile_fitted->IsOpen()){
              runfile->Close(); delete runfile; runfile = 0;
              runfile = TFile::Open(filename_nonfitted.data(), "READ");
              runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
              if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
            }
            if(read_fit){
              f0 = ReadFit(&hdQds)[0];
            }
            else{
              f0 = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, .5)[0];
              runfile_fitted->cd();
              hdQds.Write();
            }
            if(f0 > 0){
              fill_2d_graph(mpv_field_ByLEMs_view0[lem], field_it.first.first, field_it.first.second, f0);
              fill_2d_hist(h_mpv_field_ByLEMs_view0[lem], field_it.first.first, field_it.first.second, f0);
            }
            
            name_to_get = "dQds_LEM_"+to_string(lem)+"_view1" + corrected_or_not;
            if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
            if(!read_fit && !runfile_fitted->IsOpen()){
              runfile->Close(); delete runfile; runfile = 0;
              runfile = TFile::Open(filename_nonfitted.data(), "READ");
              runfile_fitted = TFile::Open(filename_fitted.data(),"RECREATE");
              if(!get_histo_in_inputfile(hdQds, runfile, name_to_get, read_fit)){return;}
            }
            
            if(read_fit){
              f1 = ReadFit(&hdQds)[0];
            }
            else{
              f1 = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, .5)[0];
            }
            if(f1 > 0){
              fill_2d_graph(mpv_field_ByLEMs_view1[lem], field_it.first.first, field_it.first.second, f1);
              fill_2d_hist(h_mpv_field_ByLEMs_view1[lem], field_it.first.first, field_it.first.second, f1);
            }
            if(f0 > 0 && f1 > 0){
              fill_2d_graph(mpv_field_ByLEMs[lem], field_it.first.first, field_it.first.second, f0+f1);
              fill_2d_hist(h_mpv_field_ByLEMs[lem], field_it.first.first, field_it.first.second, f0+f1);
            }
          }//for lem
          if(runfile_fitted->IsOpen()){
            runfile_fitted->Close();
          }
          delete runfile_fitted;
          #if verbose
          cout << "      done reading histograms and filling graphs" << endl;
          #endif
          runfile->Close();
          delete runfile;
          
          #if verbose
          cout << "      next run" << endl;
          #endif
        }//for runs
        if(files_not_found == field_it.second.size()){
          #if verbose
          cout << "    No file for fields " << field_it.first.first << " and " << field_it.first.second << " in scan " << scan_type << " " << scan_num << "." << endl;
          #endif
          empty_fields++;
        }
        #if verbose
        cout << "    next fields" << endl;
        #endif
      }//while field2run
      if(empty_fields == field2run.size()){
        #if verbose
        cout << "  No files for scan " << scan_type << " " << scan_num << "." << endl;
        #endif
        empty_scan_num++;
        continue;
      }
      
      #if verbose
      cout << "  Writing file: " << string(ofile->GetName()).data() << endl;
      #endif
      
      if(mpv_field->GetN() == 0){
        #if verbose
        cout << "  No point in MPV/Gain graph " << mpv_field->GetName() << endl;
        #endif
        delete mpv_field; delete mpv_field_view0; delete mpv_field_view1;
      }
      else{
        //write the graphs on the root file
        #if verbose
        cout << "  Writing MPV/Gain graph " << mpv_field->GetName() << endl;
        #endif
        ofile->cd();
        draw_hist_2d(h_mpv_field, scan_type, ofile, outpath);
        draw_hist_2d(h_mpv_field_view0, scan_type, ofile, outpath);
        draw_hist_2d(h_mpv_field_view1, scan_type, ofile, outpath);
      }
      for( auto lem : lems ){
        if( mpv_field_ByLEMs[lem]->GetN() == 0 ){
          #if verbose
          cout << "  No point in MPV/Gain graph for LEM " << lem << endl;
          #endif
          delete mpv_field_ByLEMs[lem]; delete mpv_field_ByLEMs_view0[lem]; delete mpv_field_ByLEMs_view1[lem];
          continue;
        }
        #if verbose
        cout << "  Writing MPV/Gain graph for LEM " << lem << endl;
        #endif
        ofile->cd();
        draw_hist_2d(h_mpv_field_ByLEMs[lem], scan_type, ofile, outpath);
        draw_hist_2d(h_mpv_field_ByLEMs_view0[lem], scan_type, ofile, outpath);
        draw_hist_2d(h_mpv_field_ByLEMs_view1[lem], scan_type, ofile, outpath);
      }//for lem
      ofile->Close();
      delete ofile;
      #if verbose
      cout << "  next scan number" << endl;
      #endif
    }//for scan num
    #if verbose
    cout << "Last scan number done" << endl;
    #endif
    if(empty_scan_num == scan_it.second.size()){
      #if verbose
      cout << "No scan for " << scan_type << endl;
      #endif
      continue;
    }
    #if verbose
    cout<< "next scan type" << endl;
    #endif
  }//for scans
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
