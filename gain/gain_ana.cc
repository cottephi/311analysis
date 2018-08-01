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


void gain_ana(string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  if(!load_run_lists()){return;}
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);
  // Set stat options
  gStyle->SetStatY(0.35);                
  // Set y-position (fraction of pad size)
  gStyle->SetStatX(0.85);                
  // Set x-position (fraction of pad size)
  gStyle->SetStatW(0.2);                
  // Set width of stat-box (fraction of pad size)
  gStyle->SetStatH(0.13);                
  // Set height of stat-box (fraction of pad size)
  
  gStyle->SetPalette(kColorPrintableOnGrey);
  bool recreate_fit_file = false;
  
  string corrected_or_not = "";
  if(cut_type.find("tg") != string::npos){corrected_or_not = "_Dx_Corrected";}
  string name_to_get = "";
  
  //****************************************************************************
  //variables definition here
  map<string, vector<int> > scans;
  
//Generated with python macro (all runs), 50% drift tolerence
//  scans["Drift"] = {10};
//  scans["Extraction"] = {2};
//  scans["Induction"] = {1,2,3,4,5};

//  scans["Amplification"] = {2,4,5};
//  scans["Amplification"] = {100};

  scans["Drift"] = {1,2,3,4,5,6,7,8,9,10,11};
//  scans["Extraction"] = {1,2,3,4,5,6,7,8};
//  scans["Amplification"] = {1,2,3,4,5,6,7};
//  scans["All"] = {};
  
//tests
//  scans["Amplification"] = {5};
//  scans["Extraction"] = {2};
//  scans["Induction"] = {4};

//Generated with python macro (june_july_2017)
//  scans["Drift"] = {1,2,3,4,5,6,7,8};
//  scans["Extraction"] = {1,2,3,4,5,6};
//  scans["Amplification"] = {1,2,3};
//  scans["Induction"] = {1,2,3,4,5};

//Generated with python macro (august_sept_2017)
//  scans["Drift"] = {1,2};
//  scans["Extraction"] = {1,2,3};
//  scans["Amplification"] = {1};

//from Elog
//  scans["Amplification"] = {1,2};
//  scans["Induction"] = {1};

  //system(string("rm "+GainAna_Output+"*.root").data());
  load_cosmics();
  
  for( auto scan_it : scans ){
    if( scan_it.first  == "All"){
      for( int i = 0; i < 63; i++ ){
        scans["All"].push_back(i+1);
      }
    }
  }
  
  for( auto scan_it : scans ){
  
    string scan_type = scan_it.first;
    string outpath = GainAna_Output + cut_type + "/" + scan_type;
    check_and_mkdir(outpath);
//    string outpath_total = outpath;
//    outpath = outpath + "gain.root"
//    TFile *ofile_total = TFile::Open(outpath.data(), "RECREATE");
//    if(!ofile_total->IsOpen()){
//      cout << "Error opening file " << outpath << endl;
//      return;
//    }
//    else{
//      #if verbose
//      cout << "  Writing file " << outpath << endl;
//      #endif
//    }

    string field_type = scan_type;
    string unit = "kV/cm";
    if(scan_type == "Drift"){
      unit = "V/cm";
    }
    if(scan_type == "All"){
      field_type = "Amplification";
    }
    TMultiGraph* mpv_field_total = new TMultiGraph();
    mpv_field_total->SetTitle(string(";"+field_type+" field ("+unit+");gain").data());
    mpv_field_total->SetName("mpv_field_total");
    map<int, TMultiGraph*> mpv_field_ByLEMs_total;
    TMultiGraph* mpv_field_total_summedviews = new TMultiGraph();
    mpv_field_total_summedviews->SetTitle(string(";"+field_type+" field ("+unit+");gain").data());
    mpv_field_total_summedviews->SetName("mpv_field_total_summedviews");
    map<int, TMultiGraph*> mpv_field_ByLEMs_total_summedviews;
    for( auto lem : lems ){
        mpv_field_ByLEMs_total[lem] = new TMultiGraph();
        mpv_field_ByLEMs_total[lem]->SetTitle(string(";"+field_type+" field ("+unit+");gain").data());
        mpv_field_ByLEMs_total[lem]->SetName(string("mpv_field_LEM_"+to_string(lem)).data());
        mpv_field_ByLEMs_total_summedviews[lem] = new TMultiGraph();
        mpv_field_ByLEMs_total_summedviews[lem]->SetTitle(string(";"+field_type+" field ("+unit+");gain").data());
        mpv_field_ByLEMs_total_summedviews[lem]->SetName(string("mpv_field_LEM_"+to_string(lem)+"_summedviews").data());
    }
    
    int empty_scan_num = 0;
    for( auto scan_num : scan_it.second ){
      
      vector<TGraphErrors*> mpv_field;
      map<int, vector<TGraphErrors*> > mpv_field_ByLEMs;
      vector<TMultiGraph*> mpv_field_AllLEMs;
      vector<int> scan_nums_for_AllLEMs;
      #if verbose 
      cout << "  Initialising graphs..." << endl;
      #endif
      bool are_graph_ok = init_graph_gain(mpv_field, mpv_field_ByLEMs, mpv_field_AllLEMs, scan_nums_for_AllLEMs, scan_type, scan_num);
      if( !are_graph_ok){
        cout << "  Error while initialising graph" << endl;
        return;
      }
      TFile *ofile = new TFile();
      outpath = GainAna_Output+cut_type+"/"+scan_type+"/"+to_string(scan_num)+"/";
      if(scan_type != "All"){
        string outfile = outpath+"gain.root";
        check_and_mkdir(outfile);
        ofile = TFile::Open(outfile.data(), "RECREATE");
      }

      //Minimum number of hits to have correct statistics
      int min_number_of_hits = 200;

      vector<int> run_list;
      vector<float> field;
      vector<string> const_fields_names;
      vector<float> const_fields_values;

      //12 drift field variation tolerance
      // Drift :
      
      if(!load_runs(scan_type, scan_num, run_list, field, const_fields_names, const_fields_values)){continue;}
      
      #if verbose
      cout << "  Doing " << scan_type << " scan number " << scan_num << endl;
      #endif
      
      //****************************************************************************
      //define an output file

      const int n_runs = run_list.size();
      
      //****************************************************************************
      //create efiled2run map

      map<float, vector<int>> field2run;

      for( int ii=0; ii < n_runs; ii++ ){
        field2run[ field[ii] ].push_back( run_list[ii] );
      }

      //****************************************************************************
      //read files and fit the distributions
      
      int empty_fields = 0;
      
      for( auto field_it : field2run ){
      
        int found_run_for_field = 0;
        string run_name;

        #if verbose
        cout << "    Processing field: " << field_it.first << endl;
        #endif

        int files_not_found = 0;
        for(auto run : field_it.second){
        
        
          if(scan_type == "Amplification" or scan_type == "All"){
            corr = get_eff(runs_and_fields[run]["Extraction"], runs_and_fields[run]["Amplification"], runs_and_fields[run]["Induction"]);
          }
          else{
            corr = 1.;
          }
        
          if(gain_corrections.find(run) != gain_corrections.end()){corr = gain_corrections[run];}

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
          vector<double> f = {-1,-1};
          
          #if verbose 
          cout << "  Reading histograms..." << endl;
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
            f = ReadFit(&hdQds, mpv_field[0], scan_type, field_it.first);
//            if(runs_and_fields[run]["Extraction"] < 1800){}
          }
          else{
            f = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, 0.5, mpv_field[0], scan_type, field_it.first);
            runfile_fitted->cd();
            hdQds.Write();
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
            f = ReadFit(&hdQds,mpv_field[1], scan_type, field_it.first);
          }
          else{
            f = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, 0.5, mpv_field[1], scan_type, field_it.first);
            runfile_fitted->cd();
            hdQds.Write();
          }
          
          if( nhits_tot < min_number_of_hits ){
            #if verbose
            cout << "    Insuficient number of hits (" << nhits_tot << ") for run " << run << ". Skipping this run. \n  next run" << endl;
            #endif
            if(runfile_fitted->IsOpen()){
              runfile_fitted->Close();
            }
            delete runfile_fitted;
            runfile->Close();
            continue;
          }
          else{
            #if verbose
            cout << "    Processing " << nhits_tot << " hits for run " << run << endl;
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
              f = ReadFit(&hdQds,mpv_field_ByLEMs[lem][0], scan_type, field_it.first);
            }
            else{
              f = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, 0.5, mpv_field_ByLEMs[lem][0], scan_type, field_it.first);
              runfile_fitted->cd();
              hdQds.Write();
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
              f = ReadFit(&hdQds,mpv_field_ByLEMs[lem][1], scan_type, field_it.first);
            }
            else{
              f = fit_dQds(&hdQds, false, min_number_of_hits, 0.05, 0.5, mpv_field_ByLEMs[lem][1], scan_type, field_it.first);
              runfile_fitted->cd();
              hdQds.Write();
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
          
          #if verbose
          cout << "    next run" << endl;
          #endif
        }//for runs
        if(files_not_found == field_it.second.size()){
          #if verbose
          cout << "  No file for field " << field_it.first << " in scan " << scan_type << " " << scan_num << "." << endl;
          #endif
          empty_fields++;
        }
        #if verbose
        cout << "    next field" << endl;
        #endif
      }//while field2run
      
      if(empty_fields == field2run.size()){
        #if verbose
        cout << "  No files for scan " << scan_type << " " << scan_num << "." << endl;
        #endif
        empty_scan_num++;
        continue;
      }
      sum_views(mpv_field, mpv_field_ByLEMs);
      fill_multigraph(mpv_field_ByLEMs, mpv_field_AllLEMs);
      
      if(scan_type != "All"){
        #if verbose
        cout << "    Writing file: " << string(ofile->GetName()).data() << endl;
        #endif
      }
      
      if(mpv_field[2]->GetN() == 0){
        #if verbose
        cout << "    No point in MPV/Gain graph for scan " << scan_type << " " << scan_num << "." << endl;
        #endif
        continue;
      }
      //write the graphs on the root file
      if(scan_type != "All"){
        #if verbose
        cout << "    Writing MPV/Gain graph for scan " << scan_type << " " << scan_num << "." << endl;
        #endif
      }
      mpv_field[0]->SetTitle(to_string(scan_num).data());
      mpv_field[1]->SetTitle(to_string(scan_num).data());
      mpv_field[2]->SetTitle(to_string(scan_num).data());
      mpv_field_total->Add(new TGraphErrors(*mpv_field[0]),"*");
      mpv_field_total->Add(new TGraphErrors(*mpv_field[1]),"*");
      mpv_field_total_summedviews->Add(new TGraphErrors(*mpv_field[2]),"*");
      mpv_field[0]->SetMarkerColor(kRed);
      mpv_field[1]->SetMarkerColor(kBlue);
      mpv_field[2]->SetMarkerColor(kGreen+1);
      
      if(scan_type != "All"){
        ofile->cd();
        if(scan_type == "Amplification"){
          //Arho~4000  Brho~200
          fit_gain(mpv_field[0],{1,4000,200});
          fit_gain(mpv_field[1],{1,4000,200});
          if(mpv_field[2]){
            fit_gain(mpv_field[2],{1,4000,200});
          }
        }
        draw_gain_graph(mpv_field[0], scan_type, scan_num, ofile, outpath);
        draw_gain_graph(mpv_field[1], scan_type, scan_num, ofile, outpath);
        if(mpv_field[2]){
          draw_gain_graph(mpv_field[2], scan_type, scan_num, ofile, outpath);
        }
      }
//      draw_gain_graph_superposed_lems(mpv_field_AllLEMs[0], scan_type, scan_nums_for_AllLEMs, ofile, outpath);
//      draw_gain_graph_superposed_lems(mpv_field_AllLEMs[1], scan_type, scan_nums_for_AllLEMs, ofile, outpath);
//      draw_gain_graph_superposed_lems(mpv_field_AllLEMs[2], scan_type, scan_nums_for_AllLEMs, ofile, outpath);
      for( auto lem : lems ){
        if(scan_type != "All"){
          #if verbose
          cout << "    Writing MPV/Gain graph for LEM " << lem << endl;
          #endif
        }
        mpv_field_ByLEMs[lem][0]->SetTitle(to_string(scan_num).data());
        mpv_field_ByLEMs[lem][1]->SetTitle(to_string(scan_num).data());
        mpv_field_ByLEMs[lem][2]->SetTitle(to_string(scan_num).data());
        mpv_field_ByLEMs_total[lem]->Add(new TGraphErrors(*mpv_field_ByLEMs[lem][0]),"*");
        mpv_field_ByLEMs_total[lem]->Add(new TGraphErrors(*mpv_field_ByLEMs[lem][1]),"*");
        mpv_field_ByLEMs_total_summedviews[lem]->Add(new TGraphErrors(*mpv_field_ByLEMs[lem][2]),"*");
        mpv_field_ByLEMs[lem][0]->SetMarkerColor(kRed);
        mpv_field_ByLEMs[lem][1]->SetMarkerColor(kBlue);
        mpv_field_ByLEMs[lem][2]->SetMarkerColor(kGreen+1);
        
        if(scan_type != "All"){
          ofile->cd();
          if(scan_type == "Amplification"){
            //Arho~4000  Brho~200
            fit_gain(mpv_field_ByLEMs[lem][0],{1,4000,200});
            fit_gain(mpv_field_ByLEMs[lem][1],{1,4000,200});
            if(mpv_field_ByLEMs[lem][2]){
              fit_gain(mpv_field_ByLEMs[lem][2],{1,4000,200});
            }
          }
          draw_gain_graph(mpv_field_ByLEMs[lem][0], scan_type, scan_num, ofile, outpath);
          draw_gain_graph(mpv_field_ByLEMs[lem][1], scan_type, scan_num, ofile, outpath);
          if(mpv_field_ByLEMs[lem][2]){
            draw_gain_graph(mpv_field_ByLEMs[lem][2], scan_type, scan_num, ofile, outpath);
          }
        }
      }//for lem
      if(scan_type != "All"){
        ofile->Close();
        delete ofile;
      }
      
      #if verbose
      cout << "    next scan number" << endl;
      #endif
    }//for scan num
    #if verbose
    cout << "  Last scan number done" << endl;
    #endif
    if(empty_scan_num == scan_it.second.size()){
      #if verbose
      cout << "  No scan for " << scan_type << endl;
      #endif
      continue;
    }
    
//    ofile_total->cd();
//    draw_gain_multigraph(mpv_field_total, scan_type, ofile_total);
//    draw_gain_multigraph(mpv_field_total_summedviews, scan_type, ofile_total, true, outpath_total);
//    
//    ofile_total->cd();
//    for( auto lem : lems ){
//      if(mpv_field_ByLEMs_total[lem]->GetListOfGraphs()){
//        draw_gain_multigraph(mpv_field_ByLEMs_total[lem], scan_type, ofile_total);
//      }
//      if(mpv_field_ByLEMs_total_summedviews[lem]->GetListOfGraphs()){
//        draw_gain_multigraph(mpv_field_ByLEMs_total_summedviews[lem], scan_type, ofile_total, true, outpath_total);
//      }
//    }// for lem
//    ofile_total->Close();
//    delete ofile_total;
    #if verbose
    cout<< "next scan type" << endl;
    #endif
    gain_corrections.clear();
  }//for scans
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
