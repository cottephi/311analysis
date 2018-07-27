//Christoph Alt, June 2018, modified by Philippe Cotte June 2018
//christoph.alt@cern.ch
//philippe.cotte@cea.fr

// general header files:
#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <TMyFileHeader.h>

using namespace std;

//main function
void Read_tree(int file0 = 0, int file1 = 0){
  
  if(!Load_Version()){return;}
  if(file1 < file0){
    int tmp = file0;
    file0 = file1;
    file1 = tmp;
  }
  vector<track> tracks_before_cuts = {};
  vector<int> broken_files = {};//794, 658};
  gErrorIgnoreLevel = kError;
  int NFiles = 1+file1-file0;
  
  const double MeVToFc = 1.6/0.236;

  //Define some necessary constants (we don't know the size of the arrays stored in the ROOT file)

  const int NMaxNumberOfGEANTParticles = 11600; //per event
  
  const int NMaxNumberOfGEANTParticlesInTPCAV = 1954;

  const int NMaxNumberOfGEANTTrajectoryStepsForAllParticles = 230764;

  std::cout << "Looping over files " << file0 << " to " << file1 << std::endl;

  vector<track> tracks;

  for(int file=0; file<NFiles; file++){
    if ( std::find(broken_files.begin(), broken_files.end(), file+file0) != broken_files.end() ){continue;}
    //Constructing the const char* with the name (and path) of the root file
    std::stringstream stringfile;
    stringfile << file+file0;
    std::string StringRootFile( path_MC + stringfile.str() + "-G4Detsim-Parser.root");

    std::cout << std::endl;
    std::cout << "Path: " << StringRootFile << std::endl;

    const char* CharRootFile = StringRootFile.c_str();

    //Checking if subrun exists
      if(!ExistTest(StringRootFile))
      {
        std::cout << "File " << file+file0 << " doesn't exist. Skip." << std::endl;
        continue;
      }

    std::cout << "Loading file " << file+file0 << "..." << std::endl;
    //Load ROOT file and tree
    TFile *rFile = new TFile(CharRootFile, "READ");
    TTree *rTree = (TTree*)rFile->Get("analysistree/anatree"); //should always be the same
    int NEntries = (int)rTree->GetEntries();
    std::cout << "File " << file+file0 << " has " << NEntries << " events..." << std::endl;

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
//    float tMCTruth_GEANT4_TrajectoryStep_Direction_Theta[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
//    float tMCTruth_GEANT4_TrajectoryStep_Direction_Phi[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];



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
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Direction_Theta",&tMCTruth_GEANT4_TrajectoryStep_Direction_Theta);
//    rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Direction_Phi",&tMCTruth_GEANT4_TrajectoryStep_Direction_Phi);

    
    //****************
    //** Event loop **
    //****************

    for(int entry=0; entry<NEntries; entry++){
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
            if(Istep > 0){
              hit dummy_hit;
              dummy_hit.run = tRun;
              dummy_hit.subrun = tSubrun;
              dummy_hit.event = tEventNumberInRun;
              dummy_hit.view = 2;
              dummy_hit.track_id = dummy_track.id;
              dummy_hit.dq = 1000*MeVToFc*tMCTruth_GEANT4_TrajectoryStep_Energy[Istep-1] - tMCTruth_GEANT4_TrajectoryStep_Energy[Istep];
              dummy_hit.sp_x = tMCTruth_GEANT4_TrajectoryStep_Point_X[Istep];
              dummy_hit.sp_y = tMCTruth_GEANT4_TrajectoryStep_Point_Y[Istep];
              dummy_hit.sp_z = tMCTruth_GEANT4_TrajectoryStep_Point_Z[Istep];
              dummy_hit.dx = TMath::Sqrt( pow(dummy_hit.sp_x - tMCTruth_GEANT4_TrajectoryStep_Point_X[Istep-1],2) + pow(dummy_hit.sp_y - tMCTruth_GEANT4_TrajectoryStep_Point_Y[Istep-1],2) + pow(dummy_hit.sp_z - tMCTruth_GEANT4_TrajectoryStep_Point_Z[Istep-1],2) );
          dummy_track.theta = get_theta(dummy_track);
          dummy_track.phi = get_phi(dummy_track);
              dummy_hit.dqdx = dummy_hit.dq/dummy_hit.dx;
              dummy_hit.dqdx_purity_corrected = correct_dx_for_lifetime(dummy_hit.dqdx, 50-dummy_hit.sp_x);
              dummy_hit.lem = find_lem(dummy_hit.sp_y, dummy_hit.sp_z);
              
              dummy_track.hits_trk.push_back(dummy_hit);
            }//Istep == 0
          
          }//hits loop
          dummy_track.nhits = dummy_track.hits_trk.size();
          tracks.push_back(dummy_track);
        }//If in AV
        Istep+=tMCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle[particle];
      }//for particle
    } //Event loop
  } //File loop
  return;
}
