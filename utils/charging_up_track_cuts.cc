
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

void read_tree_charging_up_Feb(TChain *rTree, vector<vector<track> > &tracks_charging_up, vector<int> &TimeOfEvents){
  vector<track> tracks;
  //read the tree and store all the tree information in the respective variables
  const int NMaxHitsPerEvent=100000;
  const int NMaxClustersPerEvent=10000;
  const int NMaxTracksPerEvent=1000;
  const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
  int RunDuration = 1000; //seconds;

  //Load ROOT file and tree
  int NEntries = (int)rTree->GetEntries();

  //Define variables to store the data of the ROOT file
  //Metadata
  int tRun;
  int tSubrun;
  int tEventNumberInRun;
  int tEventTimeSeconds;
//  int  tEventTimeNanoseconds;
//  char tIsData;

  //Hit variables
  int tNumberOfHits;

  //Track variables
  short tNumberOfTracks;
  float tTrack_Length[NMaxTracksPerEventTimesNViews];

  float tTrack_StartPoint_X[NMaxTracksPerEvent];
  float tTrack_StartPoint_Y[NMaxTracksPerEvent];
  float tTrack_StartPoint_Z[NMaxTracksPerEvent];
  float tTrack_EndPoint_X[NMaxTracksPerEvent];
  float tTrack_EndPoint_Y[NMaxTracksPerEvent];
  float tTrack_EndPoint_Z[NMaxTracksPerEvent];
  
  float tTrack_StartDirection_Theta[NMaxTracksPerEvent];
  float tTrack_StartDirection_Phi[NMaxTracksPerEvent];

  float tTrack_EndDirection_Theta[NMaxTracksPerEvent];
  float tTrack_EndDirection_Phi[NMaxTracksPerEvent];

  short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2];

  //Track hit variables
  float tTrack_Hit_X[NMaxHitsPerEvent];
  float tTrack_Hit_Y[NMaxHitsPerEvent];
  float tTrack_Hit_Z[NMaxHitsPerEvent];
  float tTrack_dx_3DPosition[NMaxHitsPerEvent];
  float tTrack_dx_LocalTrackDirection[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];
  float tTrack_Hit_GoodnessOfFit[NMaxHitsPerEvent];

  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("Run",&tRun);
  rTree->SetBranchAddress("Subrun",&tSubrun);
  rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  rTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
//  //Hit variables
  rTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);

//  //Track variables
  rTree->SetBranchAddress("NumberOfTracks_pmtrack",&tNumberOfTracks);
  rTree->SetBranchAddress("Track_Length_pmtrack",&tTrack_Length);

  rTree->SetBranchAddress("Track_StartPoint_X_pmtrack", &tTrack_StartPoint_X);
  rTree->SetBranchAddress("Track_StartPoint_Y_pmtrack", &tTrack_StartPoint_Y);
  rTree->SetBranchAddress("Track_StartPoint_Z_pmtrack", &tTrack_StartPoint_Z);
  rTree->SetBranchAddress("Track_EndPoint_X_pmtrack", &tTrack_EndPoint_X);
  rTree->SetBranchAddress("Track_EndPoint_Y_pmtrack", &tTrack_EndPoint_Y);
  rTree->SetBranchAddress("Track_EndPoint_Z_pmtrack", &tTrack_EndPoint_Z);

  rTree->SetBranchAddress("Track_StartDirection_Theta_pmtrack",&tTrack_StartDirection_Theta);
  rTree->SetBranchAddress("Track_StartDirection_Phi_pmtrack",&tTrack_StartDirection_Phi);

  rTree->SetBranchAddress("Track_EndDirection_Theta_pmtrack",&tTrack_EndDirection_Theta);
  rTree->SetBranchAddress("Track_EndDirection_Phi_pmtrack",&tTrack_EndDirection_Phi);

  rTree->SetBranchAddress("Track_NumberOfHitsPerView_pmtrack",&tTrack_NumberOfHitsPerView);

//  //Track hit variables
  rTree->SetBranchAddress("Track_Hit_X_pmtrack", &tTrack_Hit_X);
  rTree->SetBranchAddress("Track_Hit_Y_pmtrack", &tTrack_Hit_Y);
  rTree->SetBranchAddress("Track_Hit_Z_pmtrack", &tTrack_Hit_Z);
  rTree->SetBranchAddress("Track_Hit_dx_3DPosition_pmtrack", &tTrack_dx_3DPosition);
  rTree->SetBranchAddress("Track_Hit_dx_LocalTrackDirection_pmtrack", &tTrack_dx_LocalTrackDirection);
  rTree->SetBranchAddress("Track_Hit_ChargeSummedADC_pmtrack", &tTrack_Hit_ChargeSummedADC);
  rTree->SetBranchAddress("Track_Hit_ChargeIntegral_pmtrack", &tTrack_Hit_ChargeIntegral);
  rTree->SetBranchAddress("Track_Hit_GoodnessOfFit_pmtrack", &tTrack_Hit_GoodnessOfFit);
  bool first = true;
  int time_first = 0;
  double Efield = -1.;
  
  for(int i=0; i<NEntries; i++){ //Event loop
    rTree->GetEntry(i);
    if(i==0){
      if(runs_and_fields.size() == 0){load_run_lists();}
      if(!load_rho_run(tRun)){return;}
      Efield = runs_and_fields[tRun]["Amplification"];
    }
    if(first){
      time_first = tEventTimeSeconds;
      first = false;
    }
    if(tEventTimeSeconds-time_first >= RunDuration){
      tracks_charging_up.push_back(tracks);
      TimeOfEvents.push_back(tEventTimeSeconds - (int)(tEventTimeSeconds-time_first)/2);
      first = true;
      time_first = 0;
      tracks.clear();
    }
    //initialize classes (not pointers for the moment)

    int a=0; //Need this counter to remember where we are in this array.
    for(int j=0; j<tNumberOfTracks; j++) //Track loop
    {
      track dummy_track;

      dummy_track.run = tRun;
      dummy_track.subrun = tSubrun;
      dummy_track.event = tEventNumberInRun;
      dummy_track.id = j;
      dummy_track.start_x = tTrack_StartPoint_X[j];
      dummy_track.start_y = tTrack_StartPoint_Y[j];
      dummy_track.start_z = tTrack_StartPoint_Z[j];

      dummy_track.end_x = tTrack_EndPoint_X[j];
      dummy_track.end_y = tTrack_EndPoint_Y[j];
      dummy_track.end_z = tTrack_EndPoint_Z[j];
      
      if(fabs(tTrack_StartDirection_Theta[j]*TMath::Pi()/180.) > TMath::Pi()/2.){
        dummy_track.start_theta = fabs(tTrack_StartDirection_Theta[j]*TMath::Pi()/180. - TMath::Pi());
      }
      else{
        dummy_track.start_theta = tTrack_StartDirection_Theta[j]*TMath::Pi()/180.;
      }
      dummy_track.start_phi = tTrack_StartDirection_Phi[j]*TMath::Pi()/180.;
      
      
      if(fabs(tTrack_EndDirection_Theta[j]*TMath::Pi()/180.) > TMath::Pi()/2.){
        dummy_track.end_theta = fabs(tTrack_EndDirection_Theta[j]*TMath::Pi()/180. - TMath::Pi());
      }
      else{
        dummy_track.end_theta = tTrack_EndDirection_Theta[j]*TMath::Pi()/180.;
      }
      dummy_track.end_phi = tTrack_EndDirection_Phi[j]*TMath::Pi()/180.;

      dummy_track.length = tTrack_Length[j];

      dummy_track.theta = get_theta(dummy_track);
      dummy_track.phi = get_phi(dummy_track);
      

      for(int k=0; k<NUM_OF_VIEWS; k++) //View loop (for this track)
      {
        for(int l=a; l<a+tTrack_NumberOfHitsPerView[j][k]; l++) //Hit loop
        {
          hit dummy_hits;

          dummy_hits.run = tRun;
          dummy_hits.subrun = tSubrun;
          dummy_hits.event = tEventNumberInRun;
          dummy_hits.view = k;
          dummy_hits.track_id = dummy_track.id;
          dummy_hits.dq_sum = tTrack_Hit_ChargeSummedADC[l]/ADC2CHARGE[k];
          dummy_hits.dq_integral = tTrack_Hit_ChargeIntegral[l]/ADC2CHARGE[k];
          dummy_hits.dx_3D = tTrack_dx_3DPosition[l];
          dummy_hits.dx_local = tTrack_dx_LocalTrackDirection[l];
          dummy_hits.sp_x = tTrack_Hit_X[l];
          dummy_hits.sp_y = tTrack_Hit_Y[l];
          dummy_hits.sp_z = tTrack_Hit_Z[l];
          dummy_hits.purity_correction = correct_dx_for_lifetime(50-tTrack_Hit_X[l]);
          dummy_hits.lem  = find_lem(dummy_hits.sp_y, dummy_hits.sp_z);
          dummy_hits.gain_density_correction_factor = gain_correction_for_rho(get_density_for_hit(tEventTimeSeconds, tRun), Efield);

          dummy_track.hits_trk.push_back(dummy_hits);
        } //end hits
        a+=tTrack_NumberOfHitsPerView[j][k];
      } //end view

    dummy_track.nhits = dummy_track.hits_trk.size();
    tracks.push_back(dummy_track);
    }//end tracks
  }//Event loop
  return;
}

void read_tree_charging_up_June(TChain *rTree, vector<vector<track> > &tracks_charging_up, vector<int> &TimeOfEvents){
  vector<track> tracks;
  //read the tree and store all the tree information in the respective variables
  const int NMaxHitsPerEvent=100000;
  const int NMaxClustersPerEvent=10000;
  const int NMaxTracksPerEvent=1000;
  const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
  int RunDuration = 1000; //seconds;

  //Load ROOT file and tree
  int NEntries = (int)rTree->GetEntries();

  //Define variables to store the data of the ROOT file
  //Metadata
  int tRun;
  int tSubrun;
  int tEventNumberInRun;
  int tEventTimeSeconds;
//  int  tEventTimeNanoseconds;
//  char tIsData;

  //Hit variables
  int tNumberOfHits;

  //Track variables
  short tNumberOfTracks;
  float tTrack_Length_Trajectory[NMaxTracksPerEventTimesNViews];
  float tTrack_Length_StraightLine[NMaxTracksPerEventTimesNViews];

  float tTrack_StartPoint_X[NMaxTracksPerEvent];
  float tTrack_StartPoint_Y[NMaxTracksPerEvent];
  float tTrack_StartPoint_Z[NMaxTracksPerEvent];
  float tTrack_EndPoint_X[NMaxTracksPerEvent];
  float tTrack_EndPoint_Y[NMaxTracksPerEvent];
  float tTrack_EndPoint_Z[NMaxTracksPerEvent];
  
  float tTrack_StartDirection_Theta[NMaxTracksPerEvent];
  float tTrack_StartDirection_Phi[NMaxTracksPerEvent];

  float tTrack_EndDirection_Theta[NMaxTracksPerEvent];
  float tTrack_EndDirection_Phi[NMaxTracksPerEvent];

  short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2];

  //Track hit variables
  float tTrack_Hit_X[NMaxHitsPerEvent];
  float tTrack_Hit_Y[NMaxHitsPerEvent];
  float tTrack_Hit_Z[NMaxHitsPerEvent];
  float tTrack_dx_3DPosition[NMaxHitsPerEvent];
  float tTrack_ds_3DLocalTrackDirection[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];
  float tTrack_Hit_GoodnessOfFit[NMaxHitsPerEvent];

  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("Run",&tRun);
  rTree->SetBranchAddress("Subrun",&tSubrun);
  rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  rTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
//  //Hit variables
  rTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);

//  //Track variables
  rTree->SetBranchAddress("NumberOfTracks",&tNumberOfTracks);
  rTree->SetBranchAddress("Track_Length_Trajectory",&tTrack_Length_Trajectory);
  rTree->SetBranchAddress("Track_Length_StraightLine",&tTrack_Length_StraightLine);

  rTree->SetBranchAddress("Track_StartPoint_X", &tTrack_StartPoint_X);
  rTree->SetBranchAddress("Track_StartPoint_Y", &tTrack_StartPoint_Y);
  rTree->SetBranchAddress("Track_StartPoint_Z", &tTrack_StartPoint_Z);
  rTree->SetBranchAddress("Track_EndPoint_X", &tTrack_EndPoint_X);
  rTree->SetBranchAddress("Track_EndPoint_Y", &tTrack_EndPoint_Y);
  rTree->SetBranchAddress("Track_EndPoint_Z", &tTrack_EndPoint_Z);

  rTree->SetBranchAddress("Track_StartDirection_Theta",&tTrack_StartDirection_Theta);
  rTree->SetBranchAddress("Track_StartDirection_Phi",&tTrack_StartDirection_Phi);

  rTree->SetBranchAddress("Track_EndDirection_Theta",&tTrack_EndDirection_Theta);
  rTree->SetBranchAddress("Track_EndDirection_Phi",&tTrack_EndDirection_Phi);

  rTree->SetBranchAddress("Track_NumberOfHitsPerView",&tTrack_NumberOfHitsPerView);

//  //Track hit variables
  rTree->SetBranchAddress("Track_Hit_X", &tTrack_Hit_X);
  rTree->SetBranchAddress("Track_Hit_Y", &tTrack_Hit_Y);
  rTree->SetBranchAddress("Track_Hit_Z", &tTrack_Hit_Z);
  rTree->SetBranchAddress("Track_Hit_ds_3DPosition", &tTrack_dx_3DPosition);
  rTree->SetBranchAddress("Track_Hit_ds_LocalTrackDirection", &tTrack_ds_3DLocalTrackDirection);
  rTree->SetBranchAddress("Track_Hit_ChargeSummedADC", &tTrack_Hit_ChargeSummedADC);
  rTree->SetBranchAddress("Track_Hit_ChargeIntegral", &tTrack_Hit_ChargeIntegral);
  rTree->SetBranchAddress("Track_Hit_GoodnessOfFit", &tTrack_Hit_GoodnessOfFit);
  bool first = true;
  int time_first = 0;
  double Efield = -1.;
  
  for(int i=0; i<NEntries; i++){ //Event loop
    rTree->GetEntry(i);
    if(i==0){
      if(runs_and_fields.size() == 0){load_run_lists();}
      if(!load_rho_run(tRun)){return;}
      Efield = runs_and_fields[tRun]["Amplification"];
    }
    if(first){
      time_first = tEventTimeSeconds;
      first = false;
    }
    if(tEventTimeSeconds-time_first >= RunDuration){
      tracks_charging_up.push_back(tracks);
      TimeOfEvents.push_back(tEventTimeSeconds - (int)(tEventTimeSeconds-time_first)/2);
      first = true;
      time_first = 0;
      tracks.clear();
    }
    //initialize classes (not pointers for the moment)

    int a=0; //Need this counter to remember where we are in this array.
    for(int j=0; j<tNumberOfTracks; j++) //Track loop
    {
      track dummy_track;

      dummy_track.run = tRun;
      dummy_track.subrun = tSubrun;
      dummy_track.event = tEventNumberInRun;
      dummy_track.id = j;
      dummy_track.start_x = tTrack_StartPoint_X[j];
      dummy_track.start_y = tTrack_StartPoint_Y[j];
      dummy_track.start_z = tTrack_StartPoint_Z[j];

      dummy_track.end_x = tTrack_EndPoint_X[j];
      dummy_track.end_y = tTrack_EndPoint_Y[j];
      dummy_track.end_z = tTrack_EndPoint_Z[j];
      
      if(fabs(tTrack_StartDirection_Theta[j]*TMath::Pi()/180.) > TMath::Pi()/2.){
        dummy_track.start_theta = fabs(tTrack_StartDirection_Theta[j]*TMath::Pi()/180. - TMath::Pi());
      }
      else{
        dummy_track.start_theta = tTrack_StartDirection_Theta[j]*TMath::Pi()/180.;
      }
      dummy_track.start_phi = tTrack_StartDirection_Phi[j]*TMath::Pi()/180.;
      
      
      if(fabs(tTrack_EndDirection_Theta[j]*TMath::Pi()/180.) > TMath::Pi()/2.){
        dummy_track.end_theta = fabs(tTrack_EndDirection_Theta[j]*TMath::Pi()/180. - TMath::Pi());
      }
      else{
        dummy_track.end_theta = tTrack_EndDirection_Theta[j]*TMath::Pi()/180.;
      }
      dummy_track.end_phi = tTrack_EndDirection_Phi[j]*TMath::Pi()/180.;

      dummy_track.length = tTrack_Length_StraightLine[j];

      dummy_track.theta = get_theta(dummy_track);
      dummy_track.phi = get_phi(dummy_track);
      

      for(int k=0; k<NUM_OF_VIEWS; k++) //View loop (for this track)
      {
        for(int l=a; l<a+tTrack_NumberOfHitsPerView[j][k]; l++) //Hit loop
        {
          hit dummy_hits;

          dummy_hits.run = tRun;
          dummy_hits.subrun = tSubrun;
          dummy_hits.event = tEventNumberInRun;
          dummy_hits.view = k;
          dummy_hits.track_id = dummy_track.id;
          dummy_hits.dq_sum = tTrack_Hit_ChargeSummedADC[l]/ADC2CHARGE[k];
          dummy_hits.dq_integral = tTrack_Hit_ChargeIntegral[l]/ADC2CHARGE[k];
          dummy_hits.dx_3D = tTrack_dx_3DPosition[l];
          dummy_hits.dx_local = tTrack_ds_3DLocalTrackDirection[l];
          dummy_hits.sp_x = tTrack_Hit_X[l];
          dummy_hits.sp_y = tTrack_Hit_Y[l];
          dummy_hits.sp_z = tTrack_Hit_Z[l];
          dummy_hits.purity_correction = correct_dx_for_lifetime(50-tTrack_Hit_X[l]);
          dummy_hits.lem  = find_lem(dummy_hits.sp_y, dummy_hits.sp_z);
          dummy_hits.gain_density_correction_factor = gain_correction_for_rho(get_density_for_hit(tEventTimeSeconds, tRun), Efield);

          dummy_track.hits_trk.push_back(dummy_hits);
        } //end hits
        a+=tTrack_NumberOfHitsPerView[j][k];
      } //end view

    dummy_track.nhits = dummy_track.hits_trk.size();
    tracks.push_back(dummy_track);
    }//end tracks
  }//Event loop
  return;
}

//////////////////////////////// MAIN //////////////////////////////////////////

void charging_up_track_cuts(vector<int> run_list={840}, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D", bool is_batch = false){
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
    
    double drift = runs_and_fields[run]["Drift"];
    double amplification = runs_and_fields[run]["Extraction"];
    double extraction = runs_and_fields[run]["Amplification"];
    double induction = runs_and_fields[run]["Induction"];
    
    TChain chain("analysistree/anatree");
    string foldername = path+to_string(run)+"/";
    int n_subruns = glob(string(foldername+"*").data()).size();
    if( n_subruns == 0){
      cout << foldername << " not found or empty" << endl;
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
    
    vector<vector<track> > tracks_before_cuts_charging_up;
    vector<int> TimeOfEvents;
    #if verbose
    cout << "Run " << run << " has " << chain.GetEntries() << " events " << endl;
    cout << "Reading input file..." << endl;
    #endif
    if(path_311data.find("Feb") != string::npos){
      read_tree_charging_up_Feb(&chain, tracks_before_cuts_charging_up, TimeOfEvents);
    }
    else if(path_311data.find("June") != string::npos){
      read_tree_charging_up_June(&chain, tracks_before_cuts_charging_up, TimeOfEvents);
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

    int subrun = 0;
    for(int i = 0; i < tracks_before_cuts_charging_up.size(); i++){
      
      vector<TH1D> hdQds;
      map<int, vector<TH1D> > hdQds_ByLEMs;
      vector<TH1D> hdQds_Dx_Corrected;
      map<int, vector<TH1D> > hdQds_ByLEMs_Dx_Corrected;
      map<int, vector<TH1D> > hdQds_ByDx;
      map<int, map<int, vector<TH1D> > > hdQds_ByDx_ByLEMs;
      
      #if verbose 
      cout << "  Initialising histos..." << endl;
      #endif
      bool are_histo_ok = init_histo_track_cuts(hdQds, hdQds_ByLEMs, hdQds_ByDx, hdQds_ByDx_ByLEMs, hdQds_Dx_Corrected, hdQds_ByLEMs_Dx_Corrected, TimeOfEvents[i]);
      if( !are_histo_ok){
        cout << "  Error while initialising histo tracks" << endl;
        return;
      }

      //start cuts *************************************************************
      
      vector<track> tracks;
      #if verbose
      cout << "Selecting tracks using " << cut_type << " cuts..." << endl;
      #endif
      if(!select_tracks(cut_type, tracks_before_cuts_charging_up[i], tracks, hdQds, hdQds_ByLEMs, hdQds_ByDx, hdQds_ByDx_ByLEMs, hdQds_Dx_Corrected, hdQds_ByLEMs_Dx_Corrected)){continue;}
      cout << "Doing cut " << cut_type << endl;
      if( tracks.size() == 0){
        cout << "skipping subrun " << subrun << endl;
        continue;
      }
      
      string ofilepath_dQds = dQds_charging_up_Output;
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
      ofilepath_dQds = ofilepath_dQds + to_string(run) + "/";
      if(!ExistTest(ofilepath_dQds)){
        #if verbose
        cout << "Creating directory " << ofilepath_dQds.data() << endl;
        #endif
        mkdir(ofilepath_dQds.data(),0777);
      }
      ofilepath_dQds = ofilepath_dQds + to_string(subrun)+".root";
      TFile ofile_dQds(ofilepath_dQds.data(), "RECREATE");
      if(!ofile_dQds.IsOpen()){
        cout << " ERROR: TFile " << ofilepath_dQds.data() << " can't be created " << endl;
        continue;
      }
      #if verbose
      cout << " TFile " << ofilepath_dQds.data() << " has been created " << endl;
      #endif
      
      ofile_dQds.cd();
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

      #if verbose
      cout << "end with subrun " << subrun << endl;
      #endif
      subrun++;
    }//for tracks_before_cut
  
    #if verbose
    cout << "end with run " << run << endl;
    #endif
  }//for run

  #if verbose
  cout << "All done.." << endl;
  #endif

}//end macro
