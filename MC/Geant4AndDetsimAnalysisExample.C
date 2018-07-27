//Christoph Alt, June 2018
//christoph.alt@cern.ch

// general header files:
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <stdio.h>
#include <math.h>

// needed for ROOT routines:
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraphErrors.h>

inline bool ExistTest (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

//main function
void Geant4AndDetsimAnalysisExample(int NFiles, std::string MCC)
{

  int EventToDisplay = 5; //Set event you want to display

  //Define some necessary constants (we don't know the size of the arrays stored in the ROOT file)
  const int NMaxParticlesPerEvent=100000;
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

  const int NMaxNumberOfGEANTTrajectoryStepsPerParticle = 1e5;
  const int NMaxNumberOfGEANTTrajectoryStepsForAllParticles = 1e6;



  //Defining histograms

  //Raw waveforms
  TH2I *hRawWaveformView0 = new TH2I("hRawWaveformView0", "Raw waveforms in view 0",NumberOfChannelsView0,0,NumberOfChannelsView0,NumberOfTicks,0,NumberOfTicks);
  TH2I *hRawWaveformView1 = new TH2I("hRawWaveformView1", "Raw waveforms in view 1",NumberOfChannelsView1,NumberOfChannelsView0,NumberOfChannels,NumberOfTicks,0,NumberOfTicks);

  //Photons
  TH1I *hDetectedPhotons = new TH1I("hDetectedPhotons", "Number of detected photons per event",100,0,30000);
  TH1F *hDetectedPhotonsPMT = new TH1F("hDetectedPhotonsPMT", "Number of detected photons per PMT",NPMTs,0,NPMTs);
  TH1F *hDetectedPhotonsTime = new TH1F("hDetectedPhotonsTime", "Arrival time of all detected photons, in ns",0.01*EventDuration,EventStarTime,EventStarTime+EventDuration);

  //GEANT
  TH1F *hGEANTNumberOfParticles = new TH1F("hGEANTNumberOfParticles", "Number of GEANT particles",100,0,3000);
  TH1F *hGEANTStartPointY = new TH1F("hGEANTStartPointY", "Start point in y of all GEANT particles, in cm",100,-3000,3000);

  //GEANT in TPC AV
  TH1F *hGEANTINTPCAVNumberOfParticles = new TH1F("hGEANTINTPCAVNumberOfParticles", "Number of GEANT particles in TPC AV",100,0,300);
  TH1F *hGEANTINTPCAVStartMomentum = new TH1F("hGEANTINTPCAVStartMomentum", "Start momentum of all GEANT particles when entering the TPC AV, in GeV, in cm",100,0,30);

  //GEANT trajectory steps
  TH1F *hGEANTNumberOfTrajectoryStepsForAllParticles = new TH1F("hGEANTNumberOfTrajectoryStepsForAllParticles", "Number of GEANT trajectory steps for all particles per event",100,0,100000);
  TH1F *hGEANTNumberOfTrajectoryStepsPerParticle = new TH1F("hGEANTNumberOfTrajectoryStepsPerParticle", "Number of GEANT trajectory steps per particle",100,0,10000);
  TH1F *hGEANTTrjaectoryStepsMuondEds = new TH1F("hGEANTTrjaectoryStepsMuondEds", "Trajectory step dE/ds for muons, in MeV/cm, in cm",100,0,30);


  //***************
  //** File loop **
  //***************

  std::cout << "Looping over the first " << NFiles << " files..." << std::endl;

    for(int file=0; file<NFiles; file++) //Subrun loop
    {
      //Constructing the const char* with the name (and path) of the root file
      std::stringstream stringfile;
      stringfile << file;
      std::string StringRootFile( "/eos/experiment/wa105/offline/LArSoft/MC/" + MCC + "/ROOT/g4detsim/" + stringfile.str() + "-G4Detsim-Parser.root");

      std::cout << std::endl;
      std::cout << "Path: " << StringRootFile << std::endl;

      const char* CharRootFile = StringRootFile.c_str();

      //Checking if subrun exists
      if(!ExistTest(StringRootFile))
      {
        std::cout << "File " << file << " doesn't exist. Skip." << std::endl;
        continue;
      }

      std::cout << "Loading file " << file << "..." << std::endl;
      //Load ROOT file and tree
      TFile *rFile = new TFile(CharRootFile, "READ");
      TTree *rTree = (TTree*)rFile->Get("analysistree/anatree"); //should always be the same
      int NEntries = (int)rTree->GetEntries();
      std::cout << "File " << file << " has " << NEntries << " events..." << std::endl;

      //Define variables to store the data of the ROOT file
      //Metadata
      int tRun;
      int tSubrun;
      int tEventNumberInRun;
      int tEventTimeSeconds;
      int tEventTimeNanoseconds;
      char tIsData;

      //Raw waveforms
      int tRawWaveform_NumberOfChannels;
      int tRawWaveform_NumberOfTicks;
      int tRawWaveform_NumberOfTicksInAllChannels;
      int tRawWaveform_Channel[NumberOfChannels];
      short tRawWaveform_ADC[NMaxNumberOfTicksInAllChannels];

      //Photons
      int tNumberOfDetectedPhotons;
      float tDetectedPhoton_Channel[NMaxNumberOfDetectedPhotons];
      float tDetectedPhoton_Time[NMaxNumberOfDetectedPhotons];

      //GEANT
      int tGEANTNumberOfParticles;
      int tGEANTPDGCode[NMaxNumberOfGEANTParticles];
      float tGEANTStartPointY[NMaxNumberOfGEANTParticles];

      //GEANT in TPC AV
      int tGEANTInTPCAVNumberOfParticles;
      float tGEANTInTPCAVStartMomentum[NMaxNumberOfGEANTParticles];

      //GEANT trajectory steps
      int tNumberOfGEANTTrajectoryStepsForAllParticles;
      int tNumberOfGEANTTrajectoryStepsPerParticle[NMaxNumberOfGEANTParticles];
      float tGEANTTrajectoryStepPointX[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
      float tGEANTTrajectoryStepPointY[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
      float tGEANTTrajectoryStepPointZ[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
      float tGEANTTrajectoryStepPDGCode[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];
      float tGEANTTrajectoryStepEnergy[NMaxNumberOfGEANTTrajectoryStepsForAllParticles];



      //Link branches in the ROOT file to variables
      //Metadata
      rTree->SetBranchAddress("Run",&tRun);
      rTree->SetBranchAddress("Subrun",&tSubrun);
      rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
      rTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
      rTree->SetBranchAddress("EventTimeNanoseconds",&tEventTimeNanoseconds);
      rTree->SetBranchAddress("IsData",&tIsData);

      //Raw waveforms
      rTree->SetBranchAddress("RawWaveform_NumberOfChannels",&tRawWaveform_NumberOfChannels);
      rTree->SetBranchAddress("RawWaveform_NumberOfTicks",&tRawWaveform_NumberOfTicks);
      rTree->SetBranchAddress("RawWaveform_NumberOfTicksInAllChannels",&tRawWaveform_NumberOfTicksInAllChannels);
      rTree->SetBranchAddress("RawWaveform_Channel",&tRawWaveform_Channel);
      rTree->SetBranchAddress("RawWaveform_ADC",&tRawWaveform_ADC);

      //Photons
      rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfDetectedPhotons",&tNumberOfDetectedPhotons);
      rTree->SetBranchAddress("MCTruth_GEANT4_DetectedPhoton_Channel",&tDetectedPhoton_Channel);
      rTree->SetBranchAddress("MCTruth_GEANT4_DetectedPhoton_Time",&tDetectedPhoton_Time);

      //GEANT
      rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfParticles",&tGEANTNumberOfParticles);
      rTree->SetBranchAddress("MCTruth_GEANT4_PDGCode",&tGEANTPDGCode);
      rTree->SetBranchAddress("MCTruth_GEANT4_StartPoint_Y",&tGEANTStartPointY);

      //GEANT in TPC AV
      rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_NumberOfParticles",&tGEANTInTPCAVNumberOfParticles);
      rTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum",&tGEANTInTPCAVStartMomentum);

      //GEANT trajectory steps
      rTree->SetBranchAddress("MCTruth_GEANT4_TotalNumberOfTrajectoryStepsForAllParticles",&tNumberOfGEANTTrajectoryStepsForAllParticles);
      rTree->SetBranchAddress("MCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle",&tNumberOfGEANTTrajectoryStepsPerParticle);
      rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Point_X",&tGEANTTrajectoryStepPointX);
      rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Point_Y",&tGEANTTrajectoryStepPointY);
      rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Point_Z",&tGEANTTrajectoryStepPointZ);
      rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_PDGCode",&tGEANTTrajectoryStepPDGCode);
      rTree->SetBranchAddress("MCTruth_GEANT4_TrajectoryStep_Energy",&tGEANTTrajectoryStepEnergy);




      //****************
      //** Event loop **
      //****************

      for(int i=0; i<NEntries; i++) //Event loop
      {
        rTree->GetEntry(i);

	//Raw waveform event display
	if(i==EventToDisplay) 
	{
	  for(int j=0; j < tRawWaveform_NumberOfChannels; j++) //Channel loop
	  {
	    for(int k=0; k < tRawWaveform_NumberOfTicks; k++) //Tick loop
	    {
	      if(tRawWaveform_Channel[j] < NumberOfChannelsView0) //View 0
	      {
		hRawWaveformView0->SetBinContent(tRawWaveform_Channel[j]+1,k+1,tRawWaveform_ADC[j*tRawWaveform_NumberOfTicks+k]);
	      }
	      else //View 1
	      {
		hRawWaveformView1->SetBinContent(tRawWaveform_Channel[j]-NumberOfChannelsView0+1,k+1,tRawWaveform_ADC[j*tRawWaveform_NumberOfTicks+k]);
	      }
	    }
	  }
	} //Raw waveform event display

	//Photons
	hDetectedPhotons->Fill(tNumberOfDetectedPhotons);
	for(int j=0; j < tNumberOfDetectedPhotons; j++) //Photon loop
	{
	  hDetectedPhotonsPMT->Fill(tDetectedPhoton_Channel[j]);
	  hDetectedPhotonsTime->Fill(tDetectedPhoton_Time[j]);
	}

	//GEANT particles and trajectory steps
	hGEANTNumberOfParticles->Fill(tGEANTNumberOfParticles);
	hGEANTNumberOfTrajectoryStepsForAllParticles->Fill(tNumberOfGEANTTrajectoryStepsForAllParticles);

	int tsit=0;

	for(int j=0; j < tGEANTNumberOfParticles; j++) //Particle loop
	{
	  hGEANTStartPointY->Fill(tGEANTStartPointY[j]);
	  hGEANTNumberOfTrajectoryStepsPerParticle->Fill(tNumberOfGEANTTrajectoryStepsPerParticle[j]);

	  if(std::abs(tGEANTPDGCode[j]) == 13) //Select muons and antimuons
	  {
	    for(int k=0; k < tNumberOfGEANTTrajectoryStepsPerParticle[j]; k++) //Trajectory steps loop
	    {
	      if(k > 0 && std::abs(tGEANTTrajectoryStepPointX[tsit]) < 50 && std::abs(tGEANTTrajectoryStepPointY[tsit]) < 48 && tGEANTTrajectoryStepPointZ[tsit] > 0 && tGEANTTrajectoryStepPointZ[tsit] < 288 )
	      {

	        float ds = std::sqrt( std::pow(tGEANTTrajectoryStepPointX[tsit]-tGEANTTrajectoryStepPointX[tsit-1],2) + std::pow(tGEANTTrajectoryStepPointY[tsit]-tGEANTTrajectoryStepPointY[tsit-1],2) + std::pow(tGEANTTrajectoryStepPointZ[tsit]-tGEANTTrajectoryStepPointZ[tsit-1],2) );
 	        float dE = 1000*(tGEANTTrajectoryStepEnergy[tsit-1] - tGEANTTrajectoryStepEnergy[tsit]);
	        hGEANTTrjaectoryStepsMuondEds->Fill(dE/ds);

	      }
	      tsit++;
	    }
	  }
	  else tsit+=tNumberOfGEANTTrajectoryStepsPerParticle[j];
	}

	//GEANT in TPC AV particles
	hGEANTINTPCAVNumberOfParticles->Fill(tGEANTInTPCAVNumberOfParticles);
	for(int j=0; j < tGEANTInTPCAVNumberOfParticles; j++) //Particle in AV loop
	{
	  hGEANTINTPCAVStartMomentum->Fill(tGEANTInTPCAVStartMomentum[j]);
	}

      } //Event loop
    } //File loop

  //Open root file to save histograms

  std::string StringRootFileHistograms("Geant4AndDetsimAnalysisExampleHistograms_" + MCC + ".root");
  const char* CharRootFileHistograms = StringRootFileHistograms.c_str();

  TFile* plots = new TFile(CharRootFileHistograms,"RECREATE");

  hRawWaveformView0->Write();
  hRawWaveformView1->Write();

  hDetectedPhotons->Write();
  hDetectedPhotonsPMT->Write();
  hDetectedPhotonsTime->Write();

  hGEANTNumberOfParticles->Write();
  hGEANTStartPointY->Write();

  hGEANTINTPCAVNumberOfParticles->Write();
  hGEANTINTPCAVStartMomentum->Write();
 
  hGEANTNumberOfTrajectoryStepsForAllParticles->Write();
  hGEANTNumberOfTrajectoryStepsPerParticle->Write();
  hGEANTTrjaectoryStepsMuondEds->Write();

  plots->Close();

} //G4AndDetsimAnalysisExample
