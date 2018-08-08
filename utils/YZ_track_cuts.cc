
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


bool init_histo(map<pair<int,int>, vector<TH1D> > &hdQds, map<vector<int>, vector<TH1D> > &hdQds_ByDx, map<pair<int,int>, vector<TH1D> > &hdQds_Dx_Corrected){
  
  if( dx >= 100. ){
    #if verbose
    cout << "ERROR: can not have dx superior to total height of detector" << endl;
    #endif
    return false;
  }
  for( int Y = -(int)(lem_size/dy); Y < (int)(lem_size/dy); Y++ ){
    for( int Z = 0; Z < (int)(6*lem_size/dz); Z++ ){
      pair<int,int> YZ = make_pair(Y,Z);
      string sign = "";
      if(Y < 0){sign="m";}
      string histname = "dQds_view0_dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z);
      hdQds[YZ] = {};
      hdQds[YZ].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
      histname = "dQds_view1_dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z);
      hdQds[YZ].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
      histname = "dQds_view0_Dx_Corrected_dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z);
      hdQds_Dx_Corrected[YZ] = {};
      hdQds_Dx_Corrected[YZ].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
      histname = "dQds_view1_Dx_Corrected_dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z);
      hdQds_Dx_Corrected[YZ].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
  
      for( int X = 0; X < (int)(100/dx); X++ ){
        vector<int> XYZ = {X,Y,Z};
        histname = "dQds_view0_dx_" + to_string(X) + "_dy_" + sign + to_string(abs(Y)) +  "_dz_" + to_string(Z);
        hdQds_ByDx[XYZ] = {};
        hdQds_ByDx[XYZ].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
        histname = "dQds_view1_dx_" + to_string(X) + "_dy_" + sign + to_string(abs(Y)) + "_dz_" + to_string(Z);
        hdQds_ByDx[XYZ].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
      }
    }
  }
  return true;
}

void rec_track_dQds_YZ(track t, map<pair<int,int>, vector<TH1D> > &hdQds, map<vector<int>, vector<TH1D> > &hdQds_ByDx, map<pair<int,int>, vector<TH1D> > &hdQds_Dx_Corrected){
  for(auto h : t.hits_trk){
    if( !isGood_lem(h.lem)){continue;}
    if( !IsGoodChan(h) ){continue;}
    double dq;
    double ds;
    if(method_ds == "3D"){ds = h.dx_3D;}
    else{ds = h.dx_local;}
    if(method_dQ == "sum"){dq = h.dq_sum;}
    else{dq = h.dq_integral;}
    double dqds = h.gain_density_correction_factor * dq/ds;
    //cut on dqdx
    if( dqds <= dQdx_cut_min or dqds > 50 or h.sp_x < tpc_boundaries[0] or h.sp_x >= tpc_boundaries[1] or h.sp_y <= tpc_boundaries[2] or h.sp_y >= tpc_boundaries[3] or h.sp_z < tpc_boundaries[4] or h.sp_z >= tpc_boundaries[5] ){ continue; }
    
    
    pair<int,int> YZ = make_pair( (int)(h.sp_y/dy), (int)(h.sp_z/dz) );
    if(h.sp_y < 0){YZ.first -= 1;}
    vector<int> XYZ = {};
    XYZ.push_back( (int)((h.sp_x+tpc_boundaries[1])/dx) );
    XYZ.push_back( YZ.first );
    XYZ.push_back( YZ.second );
    hdQds[YZ][h.view].Fill(dqds);
    hdQds_Dx_Corrected[YZ][h.view].Fill(dqds*h.purity_correction);
    hdQds_ByDx[XYZ][h.view].Fill(dqds);
  }//end hits;
  return;
}

bool select_tracks_YZ(string cut_type, vector<track> tracks, vector<track> & mips, map<pair<int,int>, vector<TH1D> > &hdQds, map<vector<int>, vector<TH1D> > &hdQds_ByDx, map<pair<int,int>, vector<TH1D> > &hdQds_Dx_Corrected){
  //select particles crossing the detector in any direction
  int count_mip=0;
  
  if((only_throughgoing and only_throughgoing_x) or (only_throughgoing and only_throughgoing_y) or (only_throughgoing and only_throughgoing_z) or (only_throughgoing_z and only_throughgoing_x) or (only_throughgoing_y and only_throughgoing_x) or (only_throughgoing_z and only_throughgoing_y) ){
    cout << "  ERROR : can not have 2 only_throughgoings set to true at the same time" << endl;
    return false;
  }

  for(auto t : tracks){

    int minx = tpc_boundaries[0] + vol_cut[0];
    int maxx = tpc_boundaries[1] - vol_cut[1];
    int miny = tpc_boundaries[2] + vol_cut[2];
    int maxy = tpc_boundaries[3] - vol_cut[3];
    int minz = tpc_boundaries[4] + vol_cut[4];
    int maxz = tpc_boundaries[5] - vol_cut[5];

    double mag = t.length;
    
    if(cut_type == "before_cuts"){
      rec_track_dQds_YZ(t, hdQds, hdQds_ByDx, hdQds_Dx_Corrected);
      mips.push_back(t);
      count_mip++;
      continue;
    }
    
    if( mag > length_cut ) {
    
      if(highway){
        if(tracks_selected_by_highway.size() == 0){
          cout << "ERROR in select_track : please load highway output" << endl;
          return false;
        }
        int track_ID = 1000*(1+t.event) + t.id;
        if( find(tracks_selected_by_highway.begin(), tracks_selected_by_highway.end(), track_ID) == tracks_selected_by_highway.end() ){continue;}
      }
      //cut on the track angle. Avoid track parallel to a view, or parallel to drift direction. Also ignore bended tracks 
      if( theta_cut > 0 and (t.theta < theta_cut or t.theta > TMath::Pi()-theta_cut) ){continue;}
      if( phi_cut > 0 and (fabs(t.phi) - ((int)(fabs(t.phi)/(TMath::Pi()/2.)))*TMath::Pi()/2. < phi_cut or fabs(t.phi) - ((int)(fabs(t.phi)/(TMath::Pi()/2.)))*TMath::Pi()/2. > TMath::Pi()/2.-phi_cut) ){continue;}
      //remove hits with too small ds
      for( vector<hit>::iterator h = t.hits_trk.begin(); h != t.hits_trk.end(); ){
        double ds;
        if(method_ds == "3D"){ds = h->dx_3D;}
        else{ds = h->dx_local;}
        if( ds > ds_cut or ds < pitch or h->sp_x < minx or h->sp_x > maxx or h->sp_y < miny or h->sp_y > maxy or h->sp_z < minz or h->sp_z > maxz or h->GoF > GoF_cut ){
          h = t.hits_trk.erase(h);
        }
        else{
          ++h;
        }
      }
      
      if(!only_throughgoing and !only_throughgoing_x and !only_throughgoing_y and !only_throughgoing_z and t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds_YZ(t, hdQds, hdQds_ByDx, hdQds_Dx_Corrected);
        count_mip++;
      }
      else if( max(t.end_x, t.start_x) > maxx and min(t.end_x, t.start_x) < minx and !only_throughgoing_y and !only_throughgoing_z and t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds_YZ(t, hdQds, hdQds_ByDx, hdQds_Dx_Corrected);
        count_mip++;
      } //end if x
      else if( max(t.end_y, t.start_y) > maxy and min(t.end_y, t.start_y) < miny and !only_throughgoing_x and !only_throughgoing_z and t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds_YZ(t, hdQds, hdQds_ByDx, hdQds_Dx_Corrected);
        count_mip++;
      } //end if z
      else if( max(t.end_z, t.start_z) > maxz and min(t.end_z, t.start_z) < minz and !only_throughgoing_x and !only_throughgoing_y and t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds_YZ(t, hdQds, hdQds_ByDx, hdQds_Dx_Corrected);
        count_mip++;
      } //end if y
    }//end mag
  }//end for tracks
  
  #if verbose
  cout << " Selected " << count_mip << " tracks over " << tracks.size() << " tracks " <<  endl;
  #endif
  return true;
}

//////////////////////////////// MAIN //////////////////////////////////////////

void YZ_track_cuts(vector<int> run_list={840}, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D", bool is_batch = false){
  method_ds = m_ds;
  method_dQ = m_dQ;
  IsBatch = is_batch;
  if(!Load_Version(version)){return;}
  cut_type = set_cuts(cut_type);
  
  string path = path_wa105_311data;
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
    
    map<pair<int,int>, vector<TH1D> > hdQds;
    map<pair<int,int>, vector<TH1D> > hdQds_Dx_Corrected;
    map<vector<int>, vector<TH1D> > hdQds_ByDx;
    
    #if verbose 
    cout << "  Initialising histos..." << endl;
    #endif
    bool are_histo_ok = init_histo(hdQds, hdQds_ByDx, hdQds_Dx_Corrected);
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
    vector<track> * ptracks = &tracks;
    track single_track;
    #if verbose
    cout << "Run " << run << " has " << chain.GetEntries() << " events " << endl;
    cout << "Reading input file..." << endl;
    #endif
    if(path_311data.find("Feb") != string::npos){
      read_tree_Feb(&chain, tracks_before_cuts, tstart, tend);
    }
    else if(path_311data.find("June") != string::npos){
      read_tree_June(&chain, tracks_before_cuts, tstart, tend);
    }
    else if(path_311data.find("July") != string::npos){
      read_tree_June(&chain, tracks_before_cuts, tstart, tend);
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
    
    if(!select_tracks_YZ(cut_type, tracks_before_cuts,  tracks, hdQds, hdQds_ByDx, hdQds_Dx_Corrected)){continue;}
    cout << "Doing cut " << cut_type << endl;
    if( tracks.size() == 0){
      cout << "skipping run " << run << endl;
      continue;
    }
    
    string ofilepath_dQds = dQds_YZ_Output;
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
    ofilepath_dQds = ofilepath_dQds + to_string(run)+ ".root";
    TFile ofile_dQds(ofilepath_dQds.data(), "RECREATE");
    if(!ofile_dQds.IsOpen()){
      cout << " ERROR: TFile " << ofilepath_dQds << " can't be created " << endl;
      continue;
    }
    #if verbose
    cout << " TFile " << ofilepath_dQds << " has been created " << endl;
    #endif
    
    
    for( int Y = -(int)(lem_size/dy); Y < (int)(lem_size/dy); Y++ ){
      for( int Z = 0; Z < (int)(6*lem_size/dz); Z++ ){
        pair<int,int> YZ = make_pair(Y,Z);
        hdQds[YZ][0].Write();
        hdQds[YZ][1].Write();
        hdQds_Dx_Corrected[YZ][0].Write();
        hdQds_Dx_Corrected[YZ][1].Write();
        for( int X = 0; X < (int)(100/dx); X++ ){
          vector <int> XYZ = {X,Y,Z};
          hdQds_ByDx[XYZ][0].Write();
          hdQds_ByDx[XYZ][1].Write();
        }
      }
    }
    ofile_dQds.Close();
  }//for run
  #if verbose
  cout << "All done.." << endl;
  #endif
}//end macro
