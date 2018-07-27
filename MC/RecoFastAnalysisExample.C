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

// needed for ROOT routines:
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraphErrors.h>

inline bool ExistTest(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

//main function
void RecoFastAnalysisExample(int NFiles, std::string MCC){

  //Define some necessary constants (we don't know the size of the arrays stored in the ROOT file)
  const int NMaxHitsPerEvent=100000;
  const int NMaxClustersPerEvent=10000;
  const int NMaxTracksPerEvent=1000;
  const int NViews=2;
  const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NViews;

  int NChannelView0=320;
  int NChannelView1=960;

  double driftvelocity=0.16*1000;  //cm/mus

  double ADC2CHARGE = 67.; //to be confirmed
  double ADC2CHARGEView0 = 67.; //to be confirmed
  double ADC2CHARGEView1 = 67.; //to be confirmed

  //Track selection cuts
  int TrackXUpperCut=48; //cm
  int TrackXLowerCut=-48; //cm
  int TrackZLowerCut=52; //cm
  int TrackZUpperCut=248; //cm

  int dQdsBins = 100;

  int dQdsRange=5000/ADC2CHARGE;


  //Defining histograms
  TH1I *hNumberOfHitsAllSubruns = new TH1I("hNumberOfHitsAllSubruns", "Number of hits per event",0.1*10000,0,10000);
  TH1F *hGoodTracksLength = new TH1F("hGoodTracksLength", "Length of good tracks",0.25*400,0,400);
  TH1F *hGoodTracksPhi = new TH1F("hGoodTracksPhi", "Phi of good tracks",100,-180,180);
  TH1F *hGoodTracksTheta = new TH1F("hGoodTracksTheta", "Theta of good tracks",100,0,180);
  TH1F *hGoodTracksHits = new TH1F("hGoodTracksHits", "Hits in good tracks",0.25*(NChannelView0+NChannelView1),0,NChannelView0+NChannelView1);
  TH1F *hdQds = new TH1F("hdQds", "Good tracks: dQ/ds for all hits",100,0,dQdsRange);
  TH1F *hdQdsView0 = new TH1F("hdQdsView0", "Good tracks: dQ/ds for all hits in view 0",100,0,dQdsRange);
  TH1F *hdQdsView1 = new TH1F("hdQdsView1", "Good tracks: dQ/ds for all hits in view 1",100,0,dQdsRange);

  TH2F *hGoodTracksHitsView0VsView1 = new TH2F("hGoodTracksHitsView0VsView1", "Hits in good tracks",0.25*NChannelView1,0,NChannelView1,0.25*NChannelView0,0,NChannelView0);
  TH2F *hGoodTracksHitXvsY = new TH2F("hGoodTracksHitXvsY", "Hits in good tracks: x vs. y",120,-60,60,120,-60,60);
  TH2F *hGoodTracksHitXvsZ = new TH2F("hGoodTracksHitXvsz", "Hits in good tracks: x vs. z",320,-10,320,120,-60,60);
  TH2F *hGoodTracksHitYvsZ = new TH2F("hGoodTracksHitYvsZ", "Hits in good tracks: y vs. z",320,-10,320,120,-60,60);
  TH2F *hGoodTracksThetaVsPhi = new TH2F("hGoodTracksThetaVsPhi", "Good tracks: theta vs. phi",50,-180,180,50,0,180);
  TH2F *hGooddQdsVsDrift_LocalTrackDirection = new TH2F("hGooddQdsVsDrift_LocalTrackDirection", "Good tracks: dQ/ds vs. drift (view 0)",100,0,100/driftvelocity,dQdsBins,0,dQdsRange);
  TH2F *hGooddQdsVsDrift_3DPosition = new TH2F("hGooddQdsVsDrift_3DPosition", "Good tracks: dQ/ds vs. drift (view 0)",100,0,100/driftvelocity,dQdsBins,0,dQdsRange);
  TH2F *hGooddQdsVsDriftView0 = new TH2F("hGooddQdsVsDriftView0", "Good tracks: dQ/ds vs. drift (view 0)",100,0,100/driftvelocity,dQdsBins,0,dQdsRange);
  TH2F *hGooddQdsVsDriftView1 = new TH2F("hGooddQdsVsDriftView1", "Good tracks: dQ/ds vs. drift (view 1)",100,0,100/driftvelocity,dQdsBins,0,dQdsRange);

  TH3F *hGoodTracks3D = new TH3F("hGoodTracks3D","",320,-50,50,960,0,300,100,-50,50);


  //**************
  //** Run loop **
  //**************

  std::cout << "Looping over the first " << NFiles << " files..." << std::endl;

    for(int file=0; file<NFiles; file++) //Subrun loop
    {
      //Constructing the const char* with the name (and path) of the root file
      std::stringstream stringfile;
      stringfile << file;
      std::string StringRootFile( "/eos/experiment/wa105/offline/LArSoft/MC/" + MCC + "/ROOT/recofast/" + stringfile.str() + "-RecoFast-Parser.root");
      std::cout << std::endl;
      std::cout << "Path: " << StringRootFile << std::endl;

      const char* CharRootFile = StringRootFile.c_str();

      //Checking if subrun exists
      if(!ExistTest(StringRootFile))
      {
        std::cout << "Subrun " << subrun << " of run " << run << " doesn't exist. Skip." << std::endl;
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

      //Hit parameters
      int tNumberOfHits;
      float tHit_ChargeIntegral[NMaxHitsPerEvent];
      short tHit_TrackID[NMaxHitsPerEvent];
      short tHit_View[NMaxHitsPerEvent];

      //Cluster parameters
      short tNumberOfClusters;
      short tCluster_NumberOfHits[NMaxClustersPerEvent];

      //Track parameters
      short tNumberOfTracks;

      short tTrackID[NMaxTracksPerEvent];
      short tTrack_NumberOfHits[NMaxTracksPerEvent];
      float tTrack_Length_StraightLine[NMaxTracksPerEvent];
      float tTrack_StartPoint_X[NMaxHitsPerEvent];
      float tTrack_StartPoint_Y[NMaxHitsPerEvent];
      float tTrack_StartPoint_Z[NMaxHitsPerEvent];
      float tTrack_StartPoint_DistanceToBoundary[NMaxHitsPerEvent];
      float tTrack_EndPoint_X[NMaxHitsPerEvent];
      float tTrack_EndPoint_Y[NMaxHitsPerEvent];
      float tTrack_EndPoint_Z[NMaxHitsPerEvent];
      float tTrack_EndPoint_DistanceToBoundary[NMaxHitsPerEvent];
      float tTrack_StartDirection_Theta[NMaxTracksPerEvent];
      float tTrack_StartDirection_Phi[NMaxTracksPerEvent];
      float tTrack_StartDirection_X[NMaxTracksPerEvent];
      float tTrack_StartDirection_Y[NMaxTracksPerEvent];
      float tTrack_StartDirection_Z[NMaxTracksPerEvent];
      float tTrack_EndDirection_Theta[NMaxTracksPerEvent];
      float tTrack_EndDirection_Phi[NMaxTracksPerEvent];
      float tTrack_EndDirection_X[NMaxTracksPerEvent];
      float tTrack_EndDirection_Y[NMaxTracksPerEvent];
      float tTrack_EndDirection_Z[NMaxTracksPerEvent];

      float tTrack_PitchInViews[NMaxTracksPerEventTimesNViews];
      short tTrack_NumberOfHitsPerView[NMaxTracksPerEventTimesNViews];

      float tTrack_Hit_X[NMaxHitsPerEvent];
      float tTrack_Hit_Y[NMaxHitsPerEvent];
      float tTrack_Hit_Z[NMaxHitsPerEvent];
      float tTrack_Hit_ds_LocalTrackDirection[NMaxHitsPerEvent];
      float tTrack_Hit_ds_3DPosition[NMaxHitsPerEvent];
      short tTrack_Hit_View[NMaxHitsPerEvent];
      float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
      float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];

      //Link branches in the ROOT file with parameters
      //Metadata
      rTree->SetBranchAddress("Run",&tRun);
      rTree->SetBranchAddress("Subrun",&tSubrun);
      rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
      rTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
      rTree->SetBranchAddress("EventTimeNanoseconds",&tEventTimeNanoseconds);
      rTree->SetBranchAddress("IsData",&tIsData);

      //Hit parameters
      rTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);
      rTree->SetBranchAddress("Hit_ChargeIntegral",&tHit_ChargeIntegral);
      rTree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);
      rTree->SetBranchAddress("Hit_View",&tHit_View);
      rTree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);


      //Cluster parameters
      rTree->SetBranchAddress("NumberOfClusters",&tNumberOfClusters);
      rTree->SetBranchAddress("Cluster_NumberOfHits",&tCluster_NumberOfHits);

      //Track parameters
      rTree->SetBranchAddress("NumberOfTracks",&tNumberOfTracks);
      rTree->SetBranchAddress("TrackID",&tTrackID);
      rTree->SetBranchAddress("Track_NumberOfHits",&tTrack_NumberOfHits);
      rTree->SetBranchAddress("Track_Length_StraightLine",&tTrack_Length_StraightLine);

      rTree->SetBranchAddress("Track_StartPoint_X",&tTrack_StartPoint_X);
      rTree->SetBranchAddress("Track_StartPoint_Y",&tTrack_StartPoint_Y);
      rTree->SetBranchAddress("Track_StartPoint_Z",&tTrack_StartPoint_Z);
      rTree->SetBranchAddress("Track_StartPoint_DistanceToBoundary",&tTrack_StartPoint_DistanceToBoundary);
      rTree->SetBranchAddress("Track_EndPoint_X",&tTrack_EndPoint_X);
      rTree->SetBranchAddress("Track_EndPoint_Y",&tTrack_EndPoint_Y);
      rTree->SetBranchAddress("Track_EndPoint_Z",&tTrack_EndPoint_Z);
      rTree->SetBranchAddress("Track_EndPoint_DistanceToBoundary",&tTrack_EndPoint_DistanceToBoundary);

      rTree->SetBranchAddress("Track_StartDirection_Theta",&tTrack_StartDirection_Theta);
      rTree->SetBranchAddress("Track_StartDirection_Phi",&tTrack_StartDirection_Phi);
      rTree->SetBranchAddress("Track_StartDirection_X",&tTrack_StartDirection_X);
      rTree->SetBranchAddress("Track_StartDirection_Y",&tTrack_StartDirection_Y);
      rTree->SetBranchAddress("Track_StartDirection_Z",&tTrack_StartDirection_Z);
      rTree->SetBranchAddress("Track_EndDirection_Theta",&tTrack_EndDirection_Theta);
      rTree->SetBranchAddress("Track_EndDirection_Phi",&tTrack_EndDirection_Phi);
      rTree->SetBranchAddress("Track_EndDirection_X",&tTrack_EndDirection_X);
      rTree->SetBranchAddress("Track_EndDirection_Y",&tTrack_EndDirection_Y);
      rTree->SetBranchAddress("Track_EndDirection_Z",&tTrack_EndDirection_Z);


      rTree->SetBranchAddress("Track_PitchInViews",&tTrack_PitchInViews);
      rTree->SetBranchAddress("Track_NumberOfHitsPerView",&tTrack_NumberOfHitsPerView);

      rTree->SetBranchAddress("Track_Hit_X",&tTrack_Hit_X);
      rTree->SetBranchAddress("Track_Hit_Y",&tTrack_Hit_Y);
      rTree->SetBranchAddress("Track_Hit_Z",&tTrack_Hit_Z);
      rTree->SetBranchAddress("Track_Hit_ds_LocalTrackDirection",&tTrack_Hit_ds_LocalTrackDirection);
      rTree->SetBranchAddress("Track_Hit_ds_3DPosition",&tTrack_Hit_ds_3DPosition);
      rTree->SetBranchAddress("Track_Hit_View",&tTrack_Hit_View);
      rTree->SetBranchAddress("Track_Hit_ChargeSummedADC",&tTrack_Hit_ChargeSummedADC);
      rTree->SetBranchAddress("Track_Hit_ChargeIntegral",&tTrack_Hit_ChargeIntegral);


      //****************
      //** Event loop **
      //****************

      for(int i=0; i<NEntries; i++) //Event loop
      {
        rTree->GetEntry(i);

        hNumberOfHitsAllSubruns->Fill(tNumberOfHits);

        int b=0; //Info (dQ/ds, position, etc.) for all hits of all tracks is stored in a single array. Need this counter to remember where we are in this array.



        for(int j=0; j<tNumberOfTracks; j++) //Track loop
        {
          //Track selection. Selected tracks = good tracks.
          if( std::max(tTrack_StartPoint_X[j],tTrack_EndPoint_X[j]) > TrackXUpperCut && std::min(tTrack_StartPoint_X[j],tTrack_EndPoint_X[j]) < TrackXLowerCut &&
	      std::max(tTrack_StartPoint_Z[j],tTrack_EndPoint_Z[j]) < TrackZUpperCut && std::min(tTrack_StartPoint_Z[j],tTrack_EndPoint_Z[j]) > TrackZLowerCut &&
	      tNumberOfTracks < 10 )
          {


	    int HitCountView0a = 0;
	    int HitCountView1a = 0;

            for(int k=0; k<tNumberOfHits; k++) //Hit loop
            {
	      if(tHit_View[k] == 0 && tHit_TrackID[k] == j) HitCountView0a++;
	      if(tHit_View[k] == 1 && tHit_TrackID[k] == j) HitCountView1a++;
	    }

	    int HitCountView0b = 0;
	    int HitCountView1b = 0;

	    hGoodTracksLength->Fill(tTrack_Length_StraightLine[j]);
	    hGoodTracksPhi->Fill(tTrack_StartDirection_Phi[j]);
	    hGoodTracksTheta->Fill(tTrack_StartDirection_Theta[j]);
	    hGoodTracksThetaVsPhi->Fill(tTrack_StartDirection_Phi[j],tTrack_StartDirection_Theta[j]);
	    hGoodTracksHitsView0VsView1->Fill(tTrack_NumberOfHitsPerView[2*j+1],tTrack_NumberOfHitsPerView[2*j]);
	    hGoodTracksHits->Fill(tTrack_NumberOfHitsPerView[2*j+1]+tTrack_NumberOfHitsPerView[2*j]);

            for(int l=b; l < (b+tTrack_NumberOfHits[j]); l++) //Hit loop (for this track)
            {
	      hGoodTracks3D->Fill(tTrack_Hit_Y[l],tTrack_Hit_Z[l],tTrack_Hit_X[l]);
	      hGoodTracksHitXvsY->Fill(tTrack_Hit_Y[l],tTrack_Hit_X[l]);
	      hGoodTracksHitXvsZ->Fill(tTrack_Hit_Z[l],tTrack_Hit_X[l]);
	      hGoodTracksHitYvsZ->Fill(tTrack_Hit_Z[l],tTrack_Hit_Y[l]);

	      if(tTrack_Hit_View[l] == 0) //view 0
	      {
		HitCountView0b++;
		hdQds->Fill(tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView0));
		hGooddQdsVsDrift_LocalTrackDirection->Fill(-1*(tTrack_Hit_X[l]-50)/driftvelocity, tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_LocalTrackDirection[l]*ADC2CHARGEView0));
		hGooddQdsVsDrift_3DPosition->Fill(-1*(tTrack_Hit_X[l]-50)/driftvelocity, tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView0));
		hdQdsView0->Fill( tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView0));
		hGooddQdsVsDriftView0->Fill(-1*(tTrack_Hit_X[l]-50)/driftvelocity, tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView0));
	      }
	      if(tTrack_Hit_View[l] == 1) //view 1
	      {
	        HitCountView1b++;
		hdQds->Fill(tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView1));
		hGooddQdsVsDrift_LocalTrackDirection->Fill(-1*(tTrack_Hit_X[l]-50)/driftvelocity, tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_LocalTrackDirection[l]*ADC2CHARGEView1));
		hGooddQdsVsDrift_3DPosition->Fill(-1*(tTrack_Hit_X[l]-50)/driftvelocity, tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView1));
		hdQdsView1->Fill(tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView1));
		hGooddQdsVsDriftView1->Fill(-1*(tTrack_Hit_X[l]-50)/driftvelocity, tTrack_Hit_ChargeIntegral[l]/(tTrack_Hit_ds_3DPosition[l]*ADC2CHARGEView1));
	      }
	    } //Hit loop (for this track)
          } //Good tracks
	  b+=tTrack_NumberOfHits[j];
        } //Track loop
      } //Event loop
    } //File loop

  //Open root file to save histograms
  std::string StringRootFileHistograms( "RecoFastAnalysisExampleHistograms_" + MCC + ".root");
  const char* CharRootFileHistograms = StringRootFileHistograms.c_str();

  TFile* plots = new TFile(CharRootFileHistograms,"RECREATE");

  hNumberOfHitsAllSubruns->Write();
  hGoodTracks3D->Write();
  hGoodTracksLength->Write();
  hGoodTracksPhi->Write();
  hGoodTracksTheta->Write();
  hGoodTracksHits->Write();
  hGoodTracksHitsView0VsView1->Write();
  hGoodTracksHitXvsY->Write();
  hGoodTracksHitXvsZ->Write();
  hGoodTracksHitYvsZ->Write();
  hGoodTracksThetaVsPhi->Write();
  hGooddQdsVsDrift_LocalTrackDirection->Write();
  hGooddQdsVsDrift_3DPosition->Write();
  hGooddQdsVsDriftView0->Write();
  hGooddQdsVsDriftView1->Write();
  hdQds->Write();
  hdQdsView0->Write();
  hdQdsView1->Write();

  plots->Close();

} //dQds
