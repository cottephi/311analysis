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


void gain_ana_dQds(string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);
  // Set stat options

  if(!load_run_lists()){return;}
  
  //****************************************************************************
  //variables definition here
  map<string, vector<int> > scans;
  
  int min_number_of_hits = 200;
  
//Generated with python macro (all runs), 50% drift tolerence
//  scans["Drift"] = {1,2,3,4,5,6,7,8,9,10,11};
  scans["Drift"] = {2,3,8,10};
//  scans["Extraction"] = {1,2,3,4,5,6,7,8};
  scans["Extraction"] = {2,4};
  scans["Amplification"] = {2,4,5};
//  scans["Amplification"] = {5};
//  scans["Amplification"] = {1,2,3,4,5,6,7};
  scans["Induction"] = {1,2,3,4,5};
  
  load_cosmics();
  string outpath = GainAna_dQds_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type;
  check_and_mkdir(outpath);
  
  for( auto scan_it : scans ){
  
    string scan_type = scan_it.first;
    outpath = GainAna_Output + cut_type+"/"+scan_type;
    check_and_mkdir(outpath);

    string field_type = scan_type;
    string unit = "kV/cm";
    if(scan_type == "Drift"){
      unit = "V/cm";
    }
    if(scan_type == "All"){
      field_type = "Amplification";
    }
    map<double,map<int,TH1D> > mpv_field_total_view0;
    map<double,map<int,TH1D> > mpv_field_total_view1;
    map<double,map<int,TH1D> > mpv_field_total_summedviews;
    map<int, map<double,map<int,TH1D> > > mpv_field_ByLEMs_total_view0;
    map<int, map<double,map<int,TH1D> > > mpv_field_ByLEMs_total_view1;
    map<int, map<double,map<int,TH1D> > > mpv_field_ByLEMs_total_summedviews;
    
    map<double,map<int,TH1D> > mpv_field_Dx_Corrected_total_view0;
    map<double,map<int,TH1D> > mpv_field_Dx_Corrected_total_view1;
    map<double,map<int,TH1D> > mpv_field_Dx_Corrected_total_summedviews;
    map<int, map<double,map<int,TH1D> > > mpv_field_ByLEMs_Dx_Corrected_total_view0;
    map<int, map<double,map<int,TH1D> > > mpv_field_ByLEMs_Dx_Corrected_total_view1;
    map<int, map<double,map<int,TH1D> > > mpv_field_ByLEMs_Dx_Corrected_total_summedviews;
      
    int empty_scan_num = 0;
    for( auto scan_num : scan_it.second ){

      TFile *ofile = new TFile();
      outpath = GainAna_Output+cut_type+"/"+scan_type+"/"+to_string(scan_num);
      check_and_mkdir(outpath);
      string outfile = outpath+"/gain.root";
      ofile = TFile::Open(outfile.data(), "RECREATE");

      //Minimum number of hits to have correct statistics

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
      
        map<int, TH1D > mpv_field_view0;
        map<int, TH1D > mpv_field_view1;
        map<int, TH1D > mpv_field_summed_views;
        map<int, map<int,TH1D > > mpv_field_ByLEMs_view0;
        map<int, map<int,TH1D > > mpv_field_ByLEMs_view1;
        map<int, map<int,TH1D > > mpv_field_ByLEMs_summed_views;
        
        map<int, TH1D > mpv_field_Dx_Corrected_view0;
        map<int, TH1D > mpv_field_Dx_Corrected_view1;
        map<int, TH1D > mpv_field_Dx_Corrected_summed_views;
        map<int, map<int,TH1D > > mpv_field_ByLEMs_Dx_Corrected_view0;
        map<int, map<int,TH1D > > mpv_field_ByLEMs_Dx_Corrected_view1;
        map<int, map<int,TH1D > > mpv_field_ByLEMs_Dx_Corrected_summed_views;
      
        int found_run_for_field = 0;
        string run_name;

        #if verbose
        cout << "    Processing field: " << field_it.first << endl;
        #endif

        int files_not_found = 0;
        for(auto run : field_it.second){
      
          TH1D mpv_run_view0;
          TH1D mpv_run_view1;
          TH1D mpv_run_summed_views;
          map<int, TH1D > mpv_run_ByLEMs_view0;
          map<int, TH1D > mpv_run_ByLEMs_view1;
          map<int, TH1D > mpv_run_ByLEMs_summed_views;
          
          TH1D mpv_run_Dx_Corrected_view0;
          TH1D mpv_run_Dx_Corrected_view1;
          TH1D mpv_run_Dx_Corrected_summed_views;
          map<int, TH1D > mpv_run_ByLEMs_Dx_Corrected_view0;
          map<int, TH1D > mpv_run_ByLEMs_Dx_Corrected_view1;
          map<int, TH1D > mpv_run_ByLEMs_Dx_Corrected_summed_views;
          
          TH1D *hdQds;
        
          string ifilename = dQds_Output+cut_type+"/"+to_string(run)+".root";
          if(!ExistTest(ifilename)){
            #if verbose
            cout << "No file for run " << run << endl;
            #endif
            files_not_found++;
            continue;
          }
          TFile runfile(ifilename.data(),"READ");

          #if verbose 
          cout << "  Reading dQds histograms of run " << run << "..." << endl;
          #endif
          runfile.GetObject("dQds_view0", hdQds);
          hdQds->Scale(1/hdQds->GetEntries());
          mpv_run_view0 = *hdQds;
          runfile.GetObject("dQds_view1", hdQds);
          hdQds->Scale(1/hdQds->GetEntries());
          mpv_run_view1 = *hdQds;
          mpv_run_summed_views = mpv_run_view0;
          mpv_run_summed_views.Add(&mpv_run_view1);
          
          runfile.GetObject("dQds_view0_Dx_Corrected", hdQds);
          hdQds->Scale(1/hdQds->GetEntries());
          mpv_run_Dx_Corrected_view0 = *hdQds;
          runfile.GetObject("dQds_view1_Dx_Corrected", hdQds);
          hdQds->Scale(1/hdQds->GetEntries());
          mpv_run_Dx_Corrected_view1 = *hdQds;
          mpv_run_Dx_Corrected_summed_views = mpv_run_Dx_Corrected_view0;
          mpv_run_Dx_Corrected_summed_views.Add(&mpv_run_Dx_Corrected_view1);
          
          for( auto lem : lems ){
            runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view0").data(), hdQds);
            hdQds->Scale(1/hdQds->GetEntries());
            mpv_run_ByLEMs_view0[lem] = *hdQds;
            runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view1").data(), hdQds);
            hdQds->Scale(1/hdQds->GetEntries());
            mpv_run_ByLEMs_view1[lem] = *hdQds;
            mpv_run_ByLEMs_summed_views[lem] = mpv_run_ByLEMs_view0[lem];
            mpv_run_ByLEMs_summed_views[lem].Add(&mpv_run_ByLEMs_view1[lem]);
            
            runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view0_Dx_Corrected").data(), hdQds);
            hdQds->Scale(1/hdQds->GetEntries());
            mpv_run_ByLEMs_Dx_Corrected_view0[lem] = *hdQds;
            runfile.GetObject(string("dQds_LEM_"+to_string(lem)+"_view1_Dx_Corrected").data(), hdQds);
            hdQds->Scale(1/hdQds->GetEntries());
            mpv_run_ByLEMs_Dx_Corrected_view1[lem] = *hdQds;
            mpv_run_ByLEMs_Dx_Corrected_summed_views[lem] = mpv_run_ByLEMs_Dx_Corrected_view0[lem];
            mpv_run_ByLEMs_Dx_Corrected_summed_views[lem].Add(&mpv_run_ByLEMs_Dx_Corrected_view1[lem]);
          }
          #if verbose
          cout << "    done reading histograms" << endl;
          #endif
          delete hdQds;
          runfile.Close();
          
          if(mpv_run_view0.GetEntries() > min_number_of_hits){mpv_field_view0[run] = mpv_run_view0;}
          if(mpv_run_view1.GetEntries() > min_number_of_hits){mpv_field_view1[run] = mpv_run_view1;}
          if(mpv_run_summed_views.GetEntries() > min_number_of_hits){mpv_field_summed_views[run] = mpv_run_summed_views;}
          if(mpv_run_Dx_Corrected_view0.GetEntries() > min_number_of_hits){mpv_field_Dx_Corrected_view0[run] = mpv_run_Dx_Corrected_view0;}
          if(mpv_run_Dx_Corrected_view1.GetEntries() > min_number_of_hits){mpv_field_Dx_Corrected_view1[run] = mpv_run_Dx_Corrected_view1;}
          if(mpv_run_Dx_Corrected_summed_views.GetEntries() > min_number_of_hits){mpv_field_Dx_Corrected_summed_views[run] = mpv_run_Dx_Corrected_summed_views;}
          for( auto lem : lems ){
            if(mpv_run_ByLEMs_view0[lem].GetEntries() > min_number_of_hits){mpv_field_ByLEMs_view0[lem][run] = mpv_run_ByLEMs_view0[lem];}
            if(mpv_run_ByLEMs_view1[lem].GetEntries() > min_number_of_hits){mpv_field_ByLEMs_view0[lem][run] = mpv_run_ByLEMs_view1[lem];}
            if(mpv_run_ByLEMs_summed_views[lem].GetEntries() > min_number_of_hits){mpv_field_ByLEMs_view0[lem][run] = mpv_run_ByLEMs_summed_views[lem];}
            if(mpv_run_ByLEMs_Dx_Corrected_view0[lem].GetEntries() > min_number_of_hits){mpv_field_ByLEMs_Dx_Corrected_view0[lem][run] = mpv_run_ByLEMs_Dx_Corrected_view0[lem];}
            if(mpv_run_ByLEMs_Dx_Corrected_view1[lem].GetEntries() > min_number_of_hits){mpv_field_ByLEMs_Dx_Corrected_view1[lem][run] = mpv_run_ByLEMs_Dx_Corrected_view1[lem];}
            if(mpv_run_ByLEMs_Dx_Corrected_summed_views[lem].GetEntries() > min_number_of_hits){mpv_field_ByLEMs_Dx_Corrected_summed_views[lem][run] = mpv_run_ByLEMs_Dx_Corrected_summed_views[lem];}
          }
          #if verbose
          cout << "    next run" << endl;
          #endif
        }//for runs
        if(files_not_found == field_it.second.size()){
          #if verbose
          cout << "  No files for field " << field_it.first << " in scan " << scan_type << " " << scan_num << "." << endl;
          #endif
          empty_fields++;
          #if verbose
          cout << "    next field" << endl;
          #endif
          continue;
        }
        mpv_field_total_view0[field_it.first] = mpv_field_view0;
        mpv_field_total_view1[field_it.first] = mpv_field_view1;
        mpv_field_total_summedviews[field_it.first] = mpv_field_summed_views;
        mpv_field_Dx_Corrected_total_view0[field_it.first] = mpv_field_Dx_Corrected_view0;
        mpv_field_Dx_Corrected_total_view1[field_it.first] = mpv_field_Dx_Corrected_view1;
        mpv_field_Dx_Corrected_total_summedviews[field_it.first] = mpv_field_Dx_Corrected_summed_views;
        for( auto lem : lems ){
          mpv_field_ByLEMs_total_view0[lem][field_it.first] = mpv_field_ByLEMs_view0[lem];
          mpv_field_ByLEMs_total_view1[lem][field_it.first] = mpv_field_ByLEMs_view1[lem];
          mpv_field_ByLEMs_total_summedviews[lem][field_it.first] = mpv_field_ByLEMs_summed_views[lem];
          mpv_field_ByLEMs_Dx_Corrected_total_view0[lem][field_it.first] = mpv_field_ByLEMs_Dx_Corrected_view0[lem];
          mpv_field_ByLEMs_Dx_Corrected_total_view1[lem][field_it.first] = mpv_field_ByLEMs_Dx_Corrected_view1[lem];
          mpv_field_ByLEMs_Dx_Corrected_total_summedviews[lem][field_it.first] = mpv_field_ByLEMs_Dx_Corrected_summed_views[lem];
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
      
      #if verbose
      cout << "    Writing file: " << string(ofile->GetName()).data() << endl;
      #endif
      if(mpv_field_total_summedviews.size() == 0){
        #if verbose
        cout << "    No dQds histo for scan " << scan_type << " " << scan_num << "." << endl;
        #endif
      }
      else{
        //write the graphs on the root file
        #if verbose
        cout << "    Writing dQds histo for scan " << scan_type << " " << scan_num << "." << endl;
        #endif
        draw_gain_dQds(mpv_field_total_view0, scan_type, scan_num, ofile, outpath);
//        draw_gain_dQds(mpv_field_total_view1, scan_type, scan_num, ofile, outpath);
//        draw_gain_dQds(mpv_field_total_summedviews, scan_type, scan_num, ofile, outpath);
//        draw_gain_dQds(mpv_field_Dx_Corrected_total_view0, scan_type, scan_num, ofile, outpath);
//        draw_gain_dQds(mpv_field_Dx_Corrected_total_view1, scan_type, scan_num, ofile, outpath);
//        draw_gain_dQds(mpv_field_Dx_Corrected_total_summedviews, scan_type, scan_num, ofile, outpath);
      }
//      for( auto lem : lems ){
//        if( mpv_field_ByLEMs_total_summedviews[lem].size() == 0 ){
//          #if verbose
//          cout << "    No dQds histo graph for LEM " << lem << "for scan " << scan_type << " " << scan_num << "." << endl;
//          #endif
//          continue;
//        }
//        #if verbose
//        cout << "    Writing MPV/Gain graph for LEM " << lem << endl;
//        #endif
//        draw_gain_dQds(mpv_field_ByLEMs_total_view0[lem], scan_type, scan_num, ofile, outpath, lem);
//        draw_gain_dQds(mpv_field_ByLEMs_total_view1[lem], scan_type, scan_num, ofile, outpath, lem);
//        draw_gain_dQds(mpv_field_ByLEMs_total_summedviews[lem], scan_type, scan_num, ofile, outpath, lem);
//        draw_gain_dQds(mpv_field_ByLEMs_Dx_Corrected_total_view0[lem], scan_type, scan_num, ofile, outpath, lem);
//        draw_gain_dQds(mpv_field_ByLEMs_Dx_Corrected_total_view1[lem], scan_type, scan_num, ofile, outpath, lem);
//        draw_gain_dQds(mpv_field_ByLEMs_Dx_Corrected_total_summedviews[lem], scan_type, scan_num, ofile, outpath, lem);
//      }//for lem
      ofile->Close();
      delete ofile;
      
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
    #if verbose
    cout<< "next scan type" << endl;
    #endif
  }//for scans
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
