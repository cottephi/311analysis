
////////////////////////////////////////////////////////////////////////////////
//
// Produce compressed files wiht only the selected tracks (mips)
// Control histrograms are produced as well
//
////////////////////////////////////////////////////////////////////////////////

//Christoph Alt, June 2018, modified by Philippe Cotte June 2018
//christoph.alt@cern.ch
//philippe.cotte@cea.fr

// general header files:
#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>

using namespace std;

pair<int,int> get_strip(double z,double y){
  int strip0 = (int)z/pitch;
  int strip1 = (int)((y+tpc_boundaries[3])/pitch);
  return make_pair(strip0,strip1);
}
void read_tree_MC(TChain *rTree, vector<track> &tracks){

  
  const double MeVToFc = 1.6/0.236;

  //Define some necessary constants (we don't know the size of the arrays stored in the ROOT file)

  const int NMaxNumberOfGEANTParticles = 11600; //per event
  
  const int NMaxNumberOfGEANTParticlesInTPCAV = 1954;

  const int NMaxNumberOfGEANTTrajectoryStepsForAllParticles = 230764;

  

  //Define variables to store the data of the ROOT file
  //Metadata
  int tRun;
  int tSubrun;
  int tEventNumberInRun;

  //GEANT
  int tMCTruth_GEANT4_NumberOfParticles;//per event
  int tMCTruth_GEANT4_IsInTPCAV[NMaxNumberOfGEANTParticles];
  int tMCTruth_GEANT4_ParticleID[NMaxNumberOfGEANTParticles];

  //GEANT in TPC AV
//    int tMCTruth_GEANT4_InTPCAV_NumberOfParticles; //per event
  int tMCTruth_GEANT4_InTPCAV_PDGCode[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_StartPoint_X[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_StartPoint_Y[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_StartPoint_Z[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_StartEnergy[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_StartMomentum[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_StartMomentum_X[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_StartMomentum_Y[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_StartMomentum_Z[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_StartDirection_Theta[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_StartDirection_Phi[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_EndPoint_X[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_EndPoint_Y[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_EndPoint_Z[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_EndEnergy[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_EndMomentum[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_EndMomentum_X[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_EndMomentum_Y[NMaxNumberOfGEANTParticlesInTPCAV];
  int tMCTruth_GEANT4_InTPCAV_EndMomentum_Z[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_EndDirection_Theta[NMaxNumberOfGEANTParticlesInTPCAV];
//    int tMCTruth_GEANT4_InTPCAV_EndDirection_Phi[NMaxNumberOfGEANTParticlesInTPCAV];
  
  //GEANT trajectory steps
  int tMCTruth_GEANT4_TotalNumberOfTrajectoryStepsForAllParticles;//all steps of the entire event
  int tMCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle[NMaxNumberOfGEANTParticles];
//    int tMCTruth_GEANT4_TrajectoryStep_ParticleID[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];//starts at 1
  float tMCTruth_GEANT4_TrajectoryStep_Point_X[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
  float tMCTruth_GEANT4_TrajectoryStep_Point_Y[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
  float tMCTruth_GEANT4_TrajectoryStep_Point_Z[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
  int tMCTruth_GEANT4_TrajectoryStep_PDGCode[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
  float tMCTruth_GEANT4_TrajectoryStep_Energy[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
//    float tMCTruth_GEANT4_TrajectoryStep_Momentum[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
//    float tMCTruth_GEANT4_TrajectoryStep_Momentum_X[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
//    float tMCTruth_GEANT4_TrajectoryStep_Momentum_Y[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
//    float tMCTruth_GEANT4_TrajectoryStep_Momentum_Z[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
    float tMCTruth_GEANT4_TrajectoryStep_Direction_Theta[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
    float tMCTruth_GEANT4_TrajectoryStep_Direction_Phi[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];



  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("Run",&tRun);
  rTree->SetBranchAddress("Subrun",&tSubrun);
  rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);

  //GEANT
  rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfParticles",&tMCTruth_GEANT4_NumberOfParticles);
  rTree->SetBranchAddress("MCTruth_GEANT4_IsInTPCAV",&tMCTruth_GEANT4_IsInTPCAV);
  rTree->SetBranchAddress("MCTruth_GEANT4_ParticleID",&tMCTruth_GEANT4_ParticleID);

  //GEANT in TPC AV
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_NumberOfParticles",&tMCTruth_GEANT4_InTPCAV_NumberOfParticles);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_PDGCode",&tMCTruth_GEANT4_InTPCAV_PDGCode);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartPoint_X",&tMCTruth_GEANT4_InTPCAV_StartPoint_X);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartPoint_Y",&tMCTruth_GEANT4_InTPCAV_StartPoint_Y);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartPoint_Z",&tMCTruth_GEANT4_InTPCAV_StartPoint_Z);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartEnergy",&tMCTruth_GEANT4_InTPCAV_StartEnergy);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum",&tMCTruth_GEANT4_InTPCAV_StartMomentum);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum_X",&tMCTruth_GEANT4_InTPCAV_StartMomentum_X);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum_Y",&tMCTruth_GEANT4_InTPCAV_StartMomentum_Y);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum_Z",&tMCTruth_GEANT4_InTPCAV_StartMomentum_Z);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartDirection_Theta",&tMCTruth_GEANT4_InTPCAV_StartDirection_Theta);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartDirection_Phi",&tMCTruth_GEANT4_InTPCAV_StartDirection_Phi);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndPoint_X",&tMCTruth_GEANT4_InTPCAV_EndPoint_X);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndPoint_Y",&tMCTruth_GEANT4_InTPCAV_EndPoint_Y);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndPoint_Z",&tMCTruth_GEANT4_InTPCAV_EndPoint_Z);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndEnergy",&tMCTruth_GEANT4_InTPCAV_EndEnergy);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndMomentum",&tMCTruth_GEANT4_InTPCAV_EndMomentum);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndMomentum_X",&tMCTruth_GEANT4_InTPCAV_EndMomentum_X);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndMomentum_Y",&tMCTruth_GEANT4_InTPCAV_EndMomentum_Y);
  rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndMomentum_Z",&tMCTruth_GEANT4_InTPCAV_EndMomentum_Z);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndDirection_Theta",&tMCTruth_GEANT4_InTPCAV_EndDirection_Theta);
//    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndDirection_Phi",&tMCTruth_GEANT4_InTPCAV_EndDirection_Phi);

  //GEANT trajectory steps
  rTree->SetBranchAddress("MCTruth_GEANT4_TotalNumberOfTrajectoryStepsForAllParticles",&tMCTruth_GEANT4_TotalNumberOfTrajectoryStepsForAllParticles);
  rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle",&tMCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle);
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_ParticleID",&tMCTruth_GEANT4_TrajectoryStep_ParticleID);
  rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Point_X",&tMCTruth_GEANT4_TrajectoryStep_Point_X);
  rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Point_Y",&tMCTruth_GEANT4_TrajectoryStep_Point_Y);
  rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Point_Z",&tMCTruth_GEANT4_TrajectoryStep_Point_Z);
  rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_PDGCode",&tMCTruth_GEANT4_TrajectoryStep_PDGCode);
  rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Energy",&tMCTruth_GEANT4_TrajectoryStep_Energy);
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Momentum",&tMCTruth_GEANT4_TrajectoryStep_Momentum);
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Momentum_X",&tMCTruth_GEANT4_TrajectoryStep_Momentum_X);
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Momentum_Y",&tMCTruth_GEANT4_TrajectoryStep_Momentum_Y);
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Momentum_Z",&tMCTruth_GEANT4_TrajectoryStep_Momentum_Z);
    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Direction_Theta",&tMCTruth_GEANT4_TrajectoryStep_Direction_Theta);
    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Direction_Phi",&tMCTruth_GEANT4_TrajectoryStep_Direction_Phi);
  
  //****************
  //** Event loop **
  //****************
  for(int entry=0; entry<rTree->GetEntries(); entry++){
    rTree->GetEntry(entry);
    int Istep = 0;
    int Iav = 0;
    
    //*******************
    //** particle loop **
    //*******************
    
    for(int particle = 0; particle < tMCTruth_GEANT4_NumberOfParticles; particle++){
      
      if(tMCTruth_GEANT4_IsInTPCAV[particle] == 1){

        track dummy_track;
        dummy_track.run = tRun;
        dummy_track.subrun = tSubrun;
        dummy_track.event = tEventNumberInRun;
        dummy_track.id = tMCTruth_GEANT4_ParticleID[particle];
        
        dummy_track.start_x = tMCTruth_GEANT4_InTPCAV_StartPoint_X[Iav];
        dummy_track.start_y = tMCTruth_GEANT4_InTPCAV_StartPoint_Y[Iav];
        dummy_track.start_z = tMCTruth_GEANT4_InTPCAV_StartPoint_Z[Iav];
//          dummy_track.start_theta = tMCTruth_GEANT4_InTPCAV_StartDirection_Theta[Iav];
//          dummy_track.start_phi = tMCTruth_GEANT4_InTPCAV_StartDirection_Phi[Iav];
        double mom = TMath::Sqrt( pow(tMCTruth_GEANT4_InTPCAV_StartMomentum_X[Iav],2) + pow(tMCTruth_GEANT4_InTPCAV_StartMomentum_X[Iav],2) + pow(tMCTruth_GEANT4_InTPCAV_StartMomentum_X[Iav],2) );
        dummy_track.start_vx = tMCTruth_GEANT4_InTPCAV_StartMomentum_X[Iav]/mom;
        dummy_track.start_vy = tMCTruth_GEANT4_InTPCAV_StartMomentum_Y[Iav]/mom;
        dummy_track.start_vz = tMCTruth_GEANT4_InTPCAV_StartMomentum_Z[Iav]/mom;
        
        dummy_track.end_x = tMCTruth_GEANT4_InTPCAV_EndPoint_X[Iav];
        dummy_track.end_y = tMCTruth_GEANT4_InTPCAV_EndPoint_Y[Iav];
        dummy_track.end_z = tMCTruth_GEANT4_InTPCAV_EndPoint_Z[Iav];
//          dummy_track.end_theta = tMCTruth_GEANT4_InTPCAV_EndDirection_Theta[Iav];
//          dummy_track.end_phi = tMCTruth_GEANT4_InTPCAV_EndDirection_Phi[Iav];
        mom = TMath::Sqrt( pow(tMCTruth_GEANT4_InTPCAV_EndMomentum_X[Iav],2) + pow(tMCTruth_GEANT4_InTPCAV_EndMomentum_X[Iav],2) + pow(tMCTruth_GEANT4_InTPCAV_EndMomentum_X[Iav],2) );
        dummy_track.end_vx = tMCTruth_GEANT4_InTPCAV_EndMomentum_X[Iav]/mom;
        dummy_track.end_vy = tMCTruth_GEANT4_InTPCAV_EndMomentum_Y[Iav]/mom;
        dummy_track.end_vz = tMCTruth_GEANT4_InTPCAV_EndMomentum_Z[Iav]/mom;
        
        dummy_track.length = TMath::Sqrt( pow(dummy_track.end_x - dummy_track.start_x,2) + pow(dummy_track.end_y - dummy_track.start_y,2) + pow(dummy_track.end_z - dummy_track.start_z,2) );
        dummy_track.theta = get_theta(dummy_track);
        dummy_track.phi = get_phi(dummy_track);
  
        Iav++;
        
        //**************
        //** hit loop **
        //**************

        for(int step = 0; step < tMCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle[particle]; step++){
          if(step > 0){
            hit dummy_hit0;
            hit dummy_hit1;
            dummy_hit0.sp_y = tMCTruth_GEANT4_TrajectoryStep_Point_Y[Istep+step];
            dummy_hit0.sp_z = tMCTruth_GEANT4_TrajectoryStep_Point_Z[Istep+step];
            dummy_hit0.lem = find_lem(dummy_hit0.sp_y, dummy_hit0.sp_z);
            if(dummy_hit0.lem == -99999){continue;}
            
            double DS = TMath::Sqrt( pow(dummy_hit1.sp_x - tMCTruth_GEANT4_TrajectoryStep_Point_X[Istep+step-1],2) + pow(dummy_hit1.sp_y - tMCTruth_GEANT4_TrajectoryStep_Point_Y[Istep+step-1],2) + pow(dummy_hit1.sp_z - tMCTruth_GEANT4_TrajectoryStep_Point_Z[Istep+step-1],2) );
            double DS_0 = DS/fabs(TMath::Sin(tMCTruth_GEANT4_TrajectoryStep_Direction_Theta[Istep+step])*TMath::Sin(tMCTruth_GEANT4_TrajectoryStep_Direction_Phi[Istep+step]));
            double DS_1 = DS/fabs(TMath::Sin(tMCTruth_GEANT4_TrajectoryStep_Direction_Theta[Istep+step])*TMath::Cos(tMCTruth_GEANT4_TrajectoryStep_Direction_Phi[Istep+step]));
            
            vector<double> dss = get_dss_from_angles(tMCTruth_GEANT4_TrajectoryStep_Direction_Theta[Istep+step],tMCTruth_GEANT4_TrajectoryStep_Direction_Phi[Istep+step]);
            dummy_hit0.dx = dss[0];
            dummy_hit0.dx_MC = DS_0;
            dummy_hit0.run = tRun;
            dummy_hit0.subrun = tSubrun;
            dummy_hit0.event = tEventNumberInRun;
            dummy_hit0.view = 0;
            dummy_hit0.track_id = dummy_track.id;
            dummy_hit0.dq = 1000*MeVToFc*(tMCTruth_GEANT4_TrajectoryStep_Energy[Istep+step-1] - tMCTruth_GEANT4_TrajectoryStep_Energy[Istep+step])/2.;
            dummy_hit0.sp_x = tMCTruth_GEANT4_TrajectoryStep_Point_X[Istep+step];
            dummy_hit0.dqdx = dummy_hit0.dq/DS_0;
            dummy_hit0.dqdx_purity_corrected = correct_dx_for_lifetime(dummy_hit0.dqdx, 50-dummy_hit0.sp_x);
            pair<int,int> strips = get_strip(dummy_hit0.sp_z,dummy_hit0.sp_y);
            dummy_track.hits_trk.push_back(dummy_hit0);
            
            dummy_hit1.sp_y = tMCTruth_GEANT4_TrajectoryStep_Point_Y[Istep+step];
            dummy_hit1.sp_z = tMCTruth_GEANT4_TrajectoryStep_Point_Z[Istep+step];
            dummy_hit1.lem = find_lem(dummy_hit1.sp_y, dummy_hit1.sp_z);
            dummy_hit1.dx = dss[1];
            dummy_hit1.dx_MC = DS_1;
            dummy_hit1.run = tRun;
            dummy_hit1.subrun = tSubrun;
            dummy_hit1.event = tEventNumberInRun;
            dummy_hit1.view = 1;
            dummy_hit1.track_id = dummy_track.id;
            dummy_hit1.dq = 1000*MeVToFc*(tMCTruth_GEANT4_TrajectoryStep_Energy[Istep+step-1] - tMCTruth_GEANT4_TrajectoryStep_Energy[Istep+step])/2.;
            dummy_hit1.sp_x = tMCTruth_GEANT4_TrajectoryStep_Point_X[Istep+step];
            dummy_hit1.dqdx = dummy_hit1.dq/DS_1;
            dummy_hit1.dqdx_purity_corrected = correct_dx_for_lifetime(dummy_hit1.dqdx, 50-dummy_hit1.sp_x);
            dummy_track.hits_trk.push_back(dummy_hit1);
          }//Istep == 0
        
        }//hits loop
        dummy_track.nhits = dummy_track.hits_trk.size();
        tracks.push_back(dummy_track);
      }//If in AV
      Istep+=tMCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle[particle];
    }//for particle
  } //Event loop
  return;
}

bool init_histo(double dx, vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected, string name){
  
  if( dx >= 100. ){
    #if verbose
    cout << "ERROR: can not have dx superior to total height of detector" << endl;
    #endif
    return false;
  }
  string histname = "dQds_view0_" + name;
  hdQds.push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
  histname = "dQds_view1_" + name;
  hdQds.push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
  histname = "dQds_view0_Dx_Corrected_" + name;
  hdQds_Dx_Corrected.push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
  histname = "dQds_view1_Dx_Corrected_" + name;
  hdQds_Dx_Corrected.push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
  
  for( int X = 0; X < (int)(100/dx); X++ ){
    int x = dx*X;
    for( int lem = 0; lem < NUM_OF_GOOD_LEMS; lem++ ){
      if(x == 0){
        histname = "dQds_LEM_"+to_string(lems[lem])+"_view0_" + name;
        hdQds_ByLEMs[lems[lem]].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
        histname = "dQds_LEM_"+to_string(lems[lem])+"_view1_" + name;
        hdQds_ByLEMs[lems[lem]].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
        histname = "dQds_LEM_"+to_string(lems[lem])+"_view0_Dx_Corrected_" + name;
        hdQds_ByLEMs_Dx_Corrected[lems[lem]].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
        histname = "dQds_LEM_"+to_string(lems[lem])+"_view1_Dx_Corrected_" + name;
        hdQds_ByLEMs_Dx_Corrected[lems[lem]].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
      }
      histname = "dQds_dx_"+to_string(x)+"_LEM_"+to_string(lems[lem])+"_view0_" + name;
      hdQds_ByDx_ByLEMs[x][lems[lem]].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
      histname = "dQds_dx_"+to_string(x)+"_LEM_"+to_string(lems[lem])+"_view1_" + name;
      hdQds_ByDx_ByLEMs[x][lems[lem]].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
    }
    histname = "dQds_dx_"+to_string(x)+"_view0_" + name;
    hdQds_ByDx[x].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
    histname = "dQds_dx_"+to_string(x)+"_view1_" + name;
    hdQds_ByDx[x].push_back(TH1D(histname.data(), histname.data(), 100, 0, 50));
  }
  return true;
}

void rec_track_dQds(track t, vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected, double ds_cut, double minx, double maxx, double miny, double maxy, double minz, double maxz){

  for(auto h : t.hits_trk){
  
    //cut on dqdx
    if( h.dx > ds_cut || h.dx < pitch || h.sp_x < minx || h.sp_x > maxx || h.sp_y < miny || h.sp_y > maxy || h.sp_z < minz || h.sp_z > maxz || h.dqdx <= dQdx_cut_min || h.dqdx > 50 || h.sp_x < tpc_boundaries[0] || h.sp_x >= tpc_boundaries[1] || h.sp_y < tpc_boundaries[2] || h.sp_y > tpc_boundaries[3] || h.sp_z < tpc_boundaries[4] || h.sp_z > tpc_boundaries[5] ){ continue; }

    //cut on lems if necessary
    if( isGood_lem(h.lem) ){
      hdQds[h.view].Fill(h.dqdx);
      hdQds_ByLEMs[h.lem][h.view].Fill(h.dqdx);
      hdQds_Dx_Corrected[h.view].Fill(h.dqdx_purity_corrected);
      hdQds_ByLEMs_Dx_Corrected[h.lem][h.view].Fill(h.dqdx_purity_corrected);
      hdQds_ByDx[(int)((int)((h.sp_x+50)/dx)*dx)][h.view].Fill(h.dqdx);
      hdQds_ByDx_ByLEMs[(int)((int)((h.sp_x+50)/dx)*dx)][h.lem][h.view].Fill(h.dqdx);
    }
    else{continue;}
  }//end hits;
  return;
}

bool select_tracks(vector<track> tracks, vector<track> & mips, vector<TH1D> &hdQds_before_cuts, map<int, vector<TH1D> > &hdQds_ByLEMs_before_cuts, map<int, vector<TH1D> > &hdQds_ByDx_before_cuts, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs_before_cuts, vector<TH1D> &hdQds_Dx_Corrected_before_cuts, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected_before_cuts, vector<TH1D> &hdQds_selected_tracks, map<int, vector<TH1D> > &hdQds_ByLEMs_selected_tracks, map<int, vector<TH1D> > &hdQds_ByDx_selected_tracks, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs_selected_tracks, vector<TH1D> &hdQds_Dx_Corrected_selected_tracks,  map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected_selected_tracks, int length_cut = 0, double theta_cut = 0, double phi_cut = 0, double pitch_cut = 1000, double ds_cut = 500, bool only_throughgoing = false, bool only_throughgoing_x = false, bool only_throughgoing_y = false, bool only_throughgoing_z = false){
  //select particles crossing the detector in any direction
  int count_mip=0;
  
  if((only_throughgoing && only_throughgoing_x) || (only_throughgoing && only_throughgoing_y) || (only_throughgoing && only_throughgoing_z) || (only_throughgoing_z && only_throughgoing_x) || (only_throughgoing_y && only_throughgoing_x) || (only_throughgoing_z && only_throughgoing_y) ){
    cout << "  ERROR : can not have only_throughgoing and only_throughgoing_x set to true at the same time" << endl;
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
    
    rec_track_dQds(t, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, 0, 0, 0, 0, 0, 0, 0);
    
    if( mag > length_cut ) {
      //cut on the track angle. Avoid track parallel to a view, or parallel to drift direction. Also ignore bended tracks 
      if( t.theta < theta_cut or t.theta > TMath::Pi()-theta_cut ){continue;}
      if( fabs(t.phi) - ((int)(fabs(t.phi)/(TMath::Pi()/2.)))*TMath::Pi()/2. < phi_cut or fabs(t.phi) - ((int)(fabs(t.phi)/(TMath::Pi()/2.)))*TMath::Pi()/2. > TMath::Pi()/2.-phi_cut ){continue;}
      //ignore track if hits are too few compared to the length of the track
      if( mag/t.nhits > pitch_cut){continue;}
      //remove hits with too small ds
      
      if(!only_throughgoing && !only_throughgoing_x && !only_throughgoing_y && !only_throughgoing_z && t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds(t, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, ds_cut, minx, maxx, miny, maxy, minz, maxz);
        count_mip++;
      }
      else if( max(t.end_x, t.start_x) > maxx && min(t.end_x, t.start_x) < minx && !only_throughgoing_y && !only_throughgoing_z && t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds(t, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, ds_cut, minx, maxx, miny, maxy, minz, maxz);
        count_mip++;
      } //end if x
      else if( max(t.end_y, t.start_y) > maxy && min(t.end_y, t.start_y) < miny && !only_throughgoing_x && !only_throughgoing_z && t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds(t, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, ds_cut, minx, maxx, miny, maxy, minz, maxz);
        count_mip++;
      } //end if z
      else if( max(t.end_z, t.start_z) > maxz && min(t.end_z, t.start_z) < minz && !only_throughgoing_x && !only_throughgoing_y && t.hits_trk.size() != 0){
        mips.push_back(t);
        rec_track_dQds(t, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, ds_cut, minx, maxx, miny, maxy, minz, maxz);
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

void track_cuts_MC(int file0 = 0, int file1 = 0, const char *c_cut_type = "Ds", const char * redo_before_cuts = "false"){
  
  if(!Load_Version()){return;}
  
    
  if(file1 < file0){
    int tmp = file0;
    file0 = file1;
    file1 = tmp;
  }
  vector<int> broken_files = {};//794, 658};
  gErrorIgnoreLevel = kError;
  int NFiles = 1+file1-file0;
  
  

  std::cout << "Looping over files " << file0 << " to " << file1 << std::endl;
  
  string path = path_MC_input;
  string ofilepath = SelectTrack_MC_Output;
  string cut_type = c_cut_type;
  gStyle->SetOptFit(11111);
  gStyle->SetPalette(kRainBow);

  //getting the filename merging info togheter
  std::string filename = "";

  vector<TH1D> hdQds_before_cuts;
  map<int, vector<TH1D> > hdQds_ByLEMs_before_cuts;
  vector<TH1D> hdQds_Dx_Corrected_before_cuts;
  map<int, vector<TH1D> > hdQds_ByLEMs_Dx_Corrected_before_cuts;
  map<int, vector<TH1D> > hdQds_ByDx_before_cuts;
  map<int, map<int, vector<TH1D> > > hdQds_ByDx_ByLEMs_before_cuts;
  vector<TH1D> hdQds_selected_tracks;
  map<int, vector<TH1D> > hdQds_ByLEMs_selected_tracks;
  vector<TH1D> hdQds_Dx_Corrected_selected_tracks;
  map<int, vector<TH1D> > hdQds_ByLEMs_Dx_Corrected_selected_tracks;
  map<int, vector<TH1D> > hdQds_ByDx_selected_tracks;
  map<int, map<int, vector<TH1D> > > hdQds_ByDx_ByLEMs_selected_tracks;
  
  #if verbose 
  cout << "  Initialising histos..." << endl;
  #endif
  bool are_histo_ok = init_histo(dx, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, "before_cuts");
  if( !are_histo_ok){
    cout << "  Error while initialising histo before_cuts" << endl;
    return;
  }
  are_histo_ok = init_histo(dx, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, "selected_tracks");
  if( !are_histo_ok){
    cout << "  Error while initialising histo selected_tracks" << endl;
    return;
  }

  TChain chain("analysistree/anatree");

  for(int file=0; file<NFiles; file++){

    if ( std::find(broken_files.begin(), broken_files.end(), file+file0) != broken_files.end() ){continue;}
    //Constructing the const char* with the name (and path) of the root file
    filename = path + to_string(file+file0) + "-G4Detsim-Parser.root";
      
    #if verbose
    cout << "Adding file: " << filename << endl;
    #endif

    //Checking if subrun exists
    if(!ExistTest(filename)){
      std::cout << "File " << filename << " doesn't exist. Skip." << std::endl;
      continue;
    }
    chain.Add(filename.data());
  }//for files
  #if verbose
  cout << "Files added. " << endl;
  #endif

  //start cuts *************************************************************
  
  vector<track> tracks_before_cuts,  selected_tracks;
  track single_track_before_cuts, single_selected_track;
  #if verbose
  cout << "Found " << chain.GetEntries() << " events " << endl;
  cout << "Reading input file..." << endl;
  #endif
  read_tree_MC(&chain, tracks_before_cuts);
  #if verbose
  cout << "Selecting tracks using " << cut_type << " cuts among " << tracks_before_cuts.size() << " tracks." << endl;
  #endif
  double cut_theta = 0.1745;
  double cut_phi = 0.1745;
  if(cut_type == "throughgoing"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, 0, 500, 500, true)){return;}
  }
  else if(cut_type == "throughgoingX"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, 0, 500, 500, false, true)){return;}
  }
  else if(cut_type == "throughgoingY"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, 0, 500, 500, false, false, true)){return;}
  }
  else if(cut_type == "throughgoingZ"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, 0, 500, 500, false, false, false, true)){return;}
  }
  else if(cut_type == "theta"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, cut_theta, 0, 500, 500)){return;}
  }
  else if(cut_type == "phi"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, cut_phi, 500, 500)){return;}
  }
  else if(cut_type == "theta_phi"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, cut_theta, cut_phi, 500, 500)){return;}
  }
  else if(cut_type == "Ds"){
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, 0, 500, 1.)){return;}
  }
  else if(cut_type == "full"){//theta, phi, ds, throughgoing
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0.1745, 0.1745, 500, 1., true)){return;}
  }
  else if(cut_type == "Ds_throughgoing"){//theta, phi, ds, throughgoing
    if(!select_tracks(tracks_before_cuts, selected_tracks, hdQds_before_cuts, hdQds_ByLEMs_before_cuts, hdQds_ByDx_before_cuts, hdQds_ByDx_ByLEMs_before_cuts, hdQds_Dx_Corrected_before_cuts, hdQds_ByLEMs_Dx_Corrected_before_cuts, hdQds_selected_tracks, hdQds_ByLEMs_selected_tracks, hdQds_ByDx_selected_tracks, hdQds_ByDx_ByLEMs_selected_tracks, hdQds_Dx_Corrected_selected_tracks, hdQds_ByLEMs_Dx_Corrected_selected_tracks, 0, 0, 0, 500, 1., true)){return;}
  }
  else{
    cout << "ERROR: unknown cut type" << endl;
  }
  cout << "Doing cut " << cut_type << endl;
  //if( selected_tracks.size() == 0){return;}
  
  string ofilepath_dQds = dQds_MC_Output+cut_type+"/";
  if(!ExistTest(ofilepath_dQds)){
    #if verbose
    cout << "Creating directory " << string(ofilepath_dQds).data() << endl;
    #endif
    mkdir(string(ofilepath_dQds).data(),0777);
  }
  string ofilename = ofilepath_dQds +"dQds.root";
  TFile ofile_dQds(ofilename.data(), "RECREATE");
  if(!ofile_dQds.IsOpen()){
    cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
    return;
  }
  #if verbose
  cout << " TFile " << ofilename << " has been created " << endl;
  #endif
  
  hdQds_selected_tracks[0].Write();
  hdQds_selected_tracks[1].Write();
  hdQds_Dx_Corrected_selected_tracks[0].Write();
  hdQds_Dx_Corrected_selected_tracks[1].Write();
  for( int X = 0; X < (int)(100./dx); X++ ){
    int x = dx*X;
    for( int lem = 0; lem < NUM_OF_GOOD_LEMS; lem++ ){
      if(x == 0){
        hdQds_ByLEMs_selected_tracks[lems[lem]][0].Write();
        hdQds_ByLEMs_selected_tracks[lems[lem]][1].Write();
        hdQds_ByLEMs_Dx_Corrected_selected_tracks[lems[lem]][0].Write();
        hdQds_ByLEMs_Dx_Corrected_selected_tracks[lems[lem]][1].Write();
        
      }
      hdQds_ByDx_ByLEMs_selected_tracks[x][lems[lem]][0].Write();
      hdQds_ByDx_ByLEMs_selected_tracks[x][lems[lem]][1].Write();
    }
    hdQds_ByDx_selected_tracks[x][0].Write();
    hdQds_ByDx_selected_tracks[x][1].Write();
  }
  ofile_dQds.Close();
  if(string(redo_before_cuts) == "true"){
    cout << " Redoing dQds before cuts" << endl;
    string ofilepath_dQds_bc = dQds_MC_Output+"before_cuts/";
    if(!ExistTest(ofilepath_dQds_bc)){
      #if verbose
      cout << "Creating directory " << string(ofilepath_dQds_bc).data() << endl;
      #endif
      mkdir(string(ofilepath_dQds_bc).data(),0777);
    }
    string ofilename_bc = ofilepath_dQds_bc+"dQds.root";
    TFile ofile_dQds_bc(ofilename_bc.data(), "RECREATE");
    if(!ofile_dQds_bc.IsOpen()){
      cout << " ERROR: TFile " << ofilename_bc << " can't be created " << endl;
      return;
    }
    #if verbose
    cout << " TFile " << ofilename_bc << " has been created " << endl;
    #endif
  
    hdQds_before_cuts[0].Write();
    hdQds_before_cuts[1].Write();
    hdQds_Dx_Corrected_before_cuts[0].Write();
    hdQds_Dx_Corrected_before_cuts[1].Write();
    for( int X = 0; X < (int)(100./dx); X++ ){
      int x = dx*X;
      for( int lem = 0; lem < NUM_OF_GOOD_LEMS; lem++ ){
        if(x == 0){
          hdQds_ByLEMs_before_cuts[lems[lem]][0].Write();
          hdQds_ByLEMs_before_cuts[lems[lem]][1].Write();
          hdQds_ByLEMs_Dx_Corrected_before_cuts[lems[lem]][0].Write();
          hdQds_ByLEMs_Dx_Corrected_before_cuts[lems[lem]][1].Write();
        }
        hdQds_ByDx_ByLEMs_before_cuts[x][lems[lem]][0].Write();
        hdQds_ByDx_ByLEMs_before_cuts[x][lems[lem]][1].Write();
      }
      hdQds_ByDx_before_cuts[x][0].Write();
      hdQds_ByDx_before_cuts[x][1].Write();
    }
    ofile_dQds_bc.Close();
  }

  #if verbose
  cout << "Creating output files..." << endl;
  #endif
  if(!ExistTest(ofilepath+cut_type+"/")){
    #if verbose
    cout << "Creating directory " << string(ofilepath+cut_type+"/").data() << endl;
    #endif
    mkdir(string(ofilepath+cut_type+"/").data(),0777);
  }
  if(string(redo_before_cuts) == "true"){
    cout << " Redo tracks before cuts" << endl;
    if(!ExistTest(ofilepath+"before_cuts/")){
      #if verbose
      cout << "Creating directory " << string(ofilepath+"before_cuts/").data() << endl;
      #endif
      mkdir(string(ofilepath+"before_cuts/").data(),0777);
    }
    string ofilename_bc=ofilepath+"before_cuts/tracks.root";
    TFile ofile_bc(ofilename_bc.data(), "RECREATE");
    if(!ofile_bc.IsOpen()){
      cout << " ERROR: TFile " << ofilename_bc << " can't be created " << endl;
      return;
    }
    #if verbose
    cout << " TFile " << ofilename_bc << " has been created " << endl;
    #endif

    TTree t_before_cuts("t_before_cuts", "All reconstructed tracks");
    t_before_cuts.Branch("single_track_before_cuts", "single_track_before_cuts", &single_track_before_cuts);
    for( int i = 0; i < tracks_before_cuts.size(); i++){
      single_track_before_cuts = tracks_before_cuts.at(i);
      t_before_cuts.Fill();
    }//end loop on tracks
    ofile_bc.Write();
    ofile_bc.Close();
  }
  
  ofilename=ofilepath+cut_type+"/tracks.root";
  TFile ofile(ofilename.data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
    return;
  }
  #if verbose
  cout << " TFile " << ofilename << " has been created " << endl;
  #endif
  TTree t_selected_tracks("t_selected_tracks", "Selected tracks");
  t_selected_tracks.Branch("single_selected_track", "single_selected_track", &single_selected_track);
  for( int i = 0; i < selected_tracks.size() ; i++){
    single_selected_track = selected_tracks.at(i);
    t_selected_tracks.Fill();
  }
  ofile.Write();
  ofile.Close();

  #if verbose
  cout << "All done.." << endl;
  #endif

}//end macro


