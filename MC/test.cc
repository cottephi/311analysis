//Christoph Alt, June 2018, modified by Philippe Cotte June 2018
//christoph.alt@cern.ch
//philippe.cotte@cea.fr

#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <TMyFileHeader.h>

using namespace std;

void test(int file0 = 0, int file1 = 0){
  if(file1 < file0){
    int tmp = file0;
    file0 = file1;
    file1 = tmp;
  }
  vector<int> broken_files = {};//794, 658};
  gErrorIgnoreLevel = kError;
  int NFiles = 1+file1-file0;
  int EventToDisplay = 5; //Set event you want to display
  
  //Define some necessary constants (we don't know the size of the arrays stored in the ROOT file)
  const int NMaxParticlesPerEvent=100000;
  int tmp_NMaxParticlesPerEvent = 0;
  const int NumberOfChannels=1280;
  const int NumberOfChannelsView0=320;
  const int NumberOfChannelsView1=960;
  const int NumberOfTicks=1667;
  const int ADCRange=4096;
  const int NMaxNumberOfTicksInAllChannels = NumberOfChannels*NumberOfTicks;

  const int NMaxNumberOfDetectedPhotons=1e6;
  const int NPMTs = 5;
  const int EventDuration = 1900e3; //ns
  const int EventStarTime = -1200e3; //ns

  const int NMaxNumberOfGEANTParticles = 1e5;
  int tmp_NMaxNumberOfGEANTParticles = 0;
  int tmp_NMaxNumberOfGEANTPrimaries = 0;
  int tmp_NMaxNumberOfGEANTInTPCAVNumberOfParticles = 0;

  const int NMaxNumberOfGEANTTrajectoryStepsPerParticle = 1e5;
  int tmp_NMaxNumberOfGEANTTrajectoryStepsPerParticle = 0;
  const int NMaxNumberOfGEANTTrajectoryStepsForAllParticles = 1e6;
  int tmp_NMaxNumberOfGEANTTrajectoryStepsForAllParticles = 0;
  
  std::cout << "Looping over files " << file0 << " to " << file1 << std::endl;
  
  for(int file=0; file<NFiles; file++){
    if ( std::find(broken_files.begin(), broken_files.end(), file+file0) != broken_files.end() ){continue;}
    std::stringstream stringfile;
    stringfile << file+file0;
    std::string StringRootFile( "/eos/experiment/wa105/offline/LArSoft/MC/MC6/ROOT/g4detsim/" + stringfile.str() + "-G4Detsim-Parser.root");

    
//    Checking if subrun exists
    if(!ExistTest(StringRootFile)){
      std::cout << "File " << file+file0 << " doesn't exist. Skip." << std::endl;
      continue;
    }
    std::cout << "Loading file " << file+file0 << "..." << std::endl;
    //Load ROOT file and tree
    TFile *rFile = new TFile(StringRootFile.data(), "READ");
    TTree *rTree = (TTree*)rFile->Get("analysistree/anatree"); //should always be the same
    int NEntries = (int)rTree->GetEntries();
    std::cout << "File " << file+file0 << " has " << NEntries << " events..." << std::endl;
    
    int tRun;
    int tSubrun;
    int tEventNumberInRun;

//    //GEANT
    int tGEANTNumberOfParticles;
    int tGEANTNumberOfPrimaries;

//    //GEANT in TPC AV
    int tGEANTInTPCAVNumberOfParticles;

//    //GEANT trajectory steps
    int tNumberOfGEANTTrajectoryStepsForAllParticles;
    int tNumberOfGEANTTrajectoryStepsPerParticle[NMaxNumberOfGEANTParticles];
//    
    rTree->SetBranchAddress("Run",&tRun);
    rTree->SetBranchAddress("Subrun",&tSubrun);
    rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);

//    //GEANT
    rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfParticles",&tGEANTNumberOfParticles);
    rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfPrimaries",&tGEANTNumberOfPrimaries);

//    //GEANT in TPC AV
    rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_NumberOfParticles",&tGEANTInTPCAVNumberOfParticles);
    rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle",&tNumberOfGEANTTrajectoryStepsPerParticle);

//    //GEANT trajectory steps
    rTree->SetBranchAddress("MCTruth_GEANT4_TotalNumberOfTrajectoryStepsForAllParticles",&tNumberOfGEANTTrajectoryStepsForAllParticles);

    for(int i=0; i<NEntries; i++){
      rTree->GetEntry(i);
      if(tmp_NMaxParticlesPerEvent < tGEANTNumberOfParticles){tmp_NMaxParticlesPerEvent = tGEANTNumberOfParticles;}
      if(tmp_NMaxNumberOfGEANTInTPCAVNumberOfParticles < tGEANTInTPCAVNumberOfParticles){tmp_NMaxNumberOfGEANTInTPCAVNumberOfParticles = tGEANTInTPCAVNumberOfParticles;}
      if(tmp_NMaxNumberOfGEANTTrajectoryStepsForAllParticles < tNumberOfGEANTTrajectoryStepsForAllParticles){tmp_NMaxNumberOfGEANTTrajectoryStepsForAllParticles = tNumberOfGEANTTrajectoryStepsForAllParticles;}
    }//entries
  }//files
  cout << tmp_NMaxParticlesPerEvent << " " << tmp_NMaxNumberOfGEANTInTPCAVNumberOfParticles << " " << tmp_NMaxNumberOfGEANTTrajectoryStepsForAllParticles << endl;
  return;
}
