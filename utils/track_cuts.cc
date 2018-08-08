
////////////////////////////////////////////////////////////////////////////////
//
// Produce compressed files wiht only the selected tracks (mips)
// Control histrograms are produced as well
//
////////////////////////////////////////////////////////////////////////////////

//#include "/eos/user/p/pcotte/311analysis/lib/311Lib.cc"
#include <stdlib.h>
#include <stdio.h>
#include <311Lib.h>

using namespace std;


//////////////////////////////// MAIN //////////////////////////////////////////

void track_cuts(vector<int> run_list={840}, string cut_type = "Ds", string version = "July", string m_dQ = "sum", string m_ds = "local", bool save_tracks = true, bool is_batch = false){
  method_ds = m_ds;
  method_dQ = m_dQ;
  IsBatch = is_batch;
  if(!Load_Version(version)){return;}
  cut_type = set_cuts(cut_type);
  
  string path = path_wa105_311data;
  int to_read = 0;
  if(!load_run_lists()){return;}
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
  gStyle->SetOptFit(11111);
  gStyle->SetPalette(kRainBow);

  //getting the filename merging info togheter
  std::string filename = "";
  

  for(auto run : run_list){
  
    if(highway){if(!load_highway(run)){return;}}
  
    double drift = runs_and_fields[run]["Drift"];
    double amplification = runs_and_fields[run]["Extraction"];
    double extraction = runs_and_fields[run]["Amplification"];
    double induction = runs_and_fields[run]["Induction"];
    int tstart = 0;
    int tend = 0;
    
    vector<TH1D> hdQds;
    map<int, vector<TH1D> > hdQds_ByLEMs;
    vector<TH1D> hdQds_Dx_Corrected;
    map<int, vector<TH1D> > hdQds_ByLEMs_Dx_Corrected;
    map<int, vector<TH1D> > hdQds_ByDx;
    map<int, map<int, vector<TH1D> > > hdQds_ByDx_ByLEMs;
    
    #if verbose 
    cout << "  Initialising histos..." << endl;
    #endif
    bool are_histo_ok = init_histo_track_cuts(hdQds, hdQds_ByLEMs, hdQds_ByDx, hdQds_ByDx_ByLEMs, hdQds_Dx_Corrected, hdQds_ByLEMs_Dx_Corrected);
    if( !are_histo_ok){
      cout << "  Error while initialising histos" << endl;
      return;
    }
    #if verbose
    cout << "  ...done" << endl;
    #endif

    TChain chain("analysistree/anatree");
    string foldername = path+to_string(run)+"/";
    int n_subruns = glob(string(foldername+"*").data()).size();
    if( n_subruns == 0){
      #if verbose
      cout << foldername << " not found or empty" << endl;
      #endif
      continue;
    }
    #if verbose
    cout << "Processing " << n_subruns << " subruns in run " << run << endl;
    #endif

    string file = "";
    for(int i=0; i< n_subruns; i++){
      
      file=to_string(run)+"-"+to_string(i);
      if(file == "1009-21"){continue;}
      if(path_311data.find("Feb") != string::npos){
        file = file + "-Parser.root";
      }
      else if(path_311data.find("June") != string::npos){
        file = file + "-RecoFast-Parser.root";
      }
      else if(path_311data.find("July") != string::npos){
        file = file + "-RecoFast-Parser.root";
      }
      else{
        cout << "ERROR: unknown reconstruction version" << endl;
        return;
      }
      filename=foldername+file;
      if( ExistTest(filename) ){
        #if verbose
        cout << "Adding file: " << filename << endl;
        #endif
        chain.Add(filename.data());
      }
      else{
        #if verbose
        cout << "File " << filename << " not found." << endl;
        #endif
        continue;
      }
    }
    #if verbose
    cout << "Files added. " << endl;
    #endif

    //start cuts *************************************************************
    
    vector<track> tracks_before_cuts,  tracks;
    track single_track;
    #if verbose
    cout << "Run " << run << " has " << chain.GetEntries() << " events " << endl;
    cout << "Reading input file..." << endl;
    #endif
    if(path_311data.find("Feb") != string::npos){
      read_tree_Feb(&chain, tracks_before_cuts, tstart, tend, to_read);
    }
    else if(path_311data.find("June") != string::npos){
      read_tree_June(&chain, tracks_before_cuts, tstart, tend, to_read);
    }
    else if(path_311data.find("July") != string::npos){
      read_tree_June(&chain, tracks_before_cuts, tstart, tend, to_read);
    }
    else{
      cout << "ERROR: unknown reconstruction version" << endl;
      return;
    }
    if(tracks_before_cuts.size() == 0){
      #if verbose
      cout << "Empty run " << run << endl;
      #endif
      continue;
    }
    #if verbose
    cout << "Selecting tracks using " << cut_type << " cuts..." << endl;
    #endif
    if(!select_tracks(cut_type, tracks_before_cuts,  tracks, hdQds, hdQds_ByLEMs, hdQds_ByDx, hdQds_ByDx_ByLEMs, hdQds_Dx_Corrected, hdQds_ByLEMs_Dx_Corrected)){continue;}
    cout << "Doing cut " << cut_type << endl;
    if( tracks.size() == 0){
      cout << "skipping run " << run << endl;
      continue;
    }
  
    string ofilepath_dQds = dQds_Output;
    if(!ExistTest(ofilepath_dQds)){
      #if verbose
      cout << "Creating directory " << ofilepath_dQds.data() << endl;
      #endif
      mkdir(ofilepath_dQds.data(),0777);
    }
    ofilepath_dQds = ofilepath_dQds+cut_type+"/";
    if(!ExistTest(ofilepath_dQds)){
      #if verbose
      cout << "Creating directory " << ofilepath_dQds.data() << endl;
      #endif
      mkdir(ofilepath_dQds.data(),0777);
    }
    ofilepath_dQds = ofilepath_dQds + to_string(run)+".root";
    TFile ofile_dQds(ofilepath_dQds.data(), "RECREATE");
    if(!ofile_dQds.IsOpen()){
      cout << " ERROR: TFile " << ofilepath_dQds << " can't be created " << endl;
      continue;
    }
    #if verbose
    cout << " TFile " << ofilepath_dQds << " has been created " << endl;
    #endif
    
    hdQds[0].Write();
    hdQds[1].Write();
    hdQds_Dx_Corrected[0].Write();
    hdQds_Dx_Corrected[1].Write();
    for( int X = 0; X < (int)(100./dx); X++ ){
      int x = dx*X;
      for( int lem = 0; lem < NUM_OF_GOOD_LEMS; lem++ ){
        if(x == 0){
          hdQds_ByLEMs[lems[lem]][0].Write();
          hdQds_ByLEMs[lems[lem]][1].Write();
          hdQds_ByLEMs_Dx_Corrected[lems[lem]][0].Write();
          hdQds_ByLEMs_Dx_Corrected[lems[lem]][1].Write();
        }
        hdQds_ByDx_ByLEMs[x][lems[lem]][0].Write();
        hdQds_ByDx_ByLEMs[x][lems[lem]][1].Write();
      }
      hdQds_ByDx[x][0].Write();
      hdQds_ByDx[x][1].Write();
    }
    ofile_dQds.Close();
    
    if(!save_tracks){return;}

    #if verbose
    cout << "Creating tracks files..." << endl;
    #endif
    string ofilepath_tracks = SelectTrack_Output;
    if(!ExistTest(ofilepath_tracks)){
      #if verbose
      cout << "Creating directory " << ofilepath_tracks.data() << endl;
      #endif
      mkdir(ofilepath_tracks.data(),0777);
    }
    ofilepath_tracks = ofilepath_tracks+cut_type+"/";
    if(!ExistTest(ofilepath_tracks)){
      #if verbose
      cout << "Creating directory " << ofilepath_tracks.data() << endl;
      #endif
      mkdir(ofilepath_tracks.data(),0777);
    }
    ofilepath_tracks=ofilepath_tracks+"tracks_"+to_string(run)+".root";
    TFile ofile(ofilepath_tracks.data(), "RECREATE");
    if(!ofile.IsOpen()){
      cout << " ERROR: TFile " << ofilepath_tracks << " can't be created " << endl;
      continue;
    }
    #if verbose
    cout << " TFile " << ofilepath_tracks << " has been created " << endl;
    #endif
    TTree t_tracks("tracks", "tracks");
    t_tracks.Branch("track", "track", &single_track);
    for(auto tr : tracks){
      single_track = tr;
      if(dray_miti){drays_mitigation(single_track);}
      t_tracks.Fill();
    }
    ofile.Write();
    ofile.Close();

    #if verbose
    cout << "end with run " << run << endl;
    #endif
  }//for run_list

  #if verbose
  cout << "All done.." << endl;
  #endif

}//end macro
