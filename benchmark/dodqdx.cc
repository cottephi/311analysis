////////////////////////////////////////////////////////////////////////////////
//
// Macro to compute the dQds and fit it with langaus distribution in a given
// detector sensible area
// and load all the functions.
// Returns an histogram with a fitted langaus function.
//
////////////////////////////////////////////////////////////////////////////////

// general header files:
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

// needed for ROOT routines:
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>

// project libraries
#include <311Lib.h>

using namespace std;

//list of cuts
  int length_cut = 100;
  double angle_cut = 3.1415926/2;

  vector<int> vol_cut = {2, 2, 2, 2, 2, 2}; //fiductial volume cut
  vector<int> lems = { 2, 4, 5, 6, 7, 8, 9, 11}; //active lems

//fit***************************************************************************
//Fit from bin with more entries up to > low
//Low is the lower limit over the fit ending point (up is 3.0*mean)
double low[NUM_OF_VIEWS] = {1.0, 1.2};

//fit from high > to bin with largest entries
//high is the upper limit over the fit starting point (lower is 0.0)
double high[NUM_OF_VIEWS] = {1.0, 0.5};


//root variables
TChain *chain;

TH1D *hdQds[NUM_OF_VIEWS];
TH1D *hRatio[NUM_OF_VIEWS];
TH2D *hRatio2D[NUM_OF_VIEWS];
TF1 *function_gain[NUM_OF_VIEWS];
TProfile *hRatioTheta[NUM_OF_VIEWS];
TProfile2D *hGainUniformity[NUM_OF_VIEWS];


void save_to_file(string ofilename){

  TFile *dummy_file = new TFile(ofilename.c_str(), "RECREATE");

  if(dummy_file->IsOpen())
    cout << "file " << ofilename << " open successfully" << endl;

  /*
  TF1 *func[2]; string fname = "";

  //fname = "Fitfcn_"+;
  func[0]= h[0]->GetFunction(fname.c_str());
  //fname = "Fitfcn_"+;
  func[1]= h[1]->GetFunction(fname.c_str());

  string cname = "cdQds_"+to_string(run);
  TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 700, 700);
  gPad->SetTickx(1); gPad->SetTicky(1);

  h[0]->SetLineColor(kRed);
  h[1]->SetLineColor(kBlue);

  h[0]->SetLineWidth(2);
  h[1]->SetLineWidth(2);

  h[0]->SetMarkerSize(0.7);
  h[1]->SetMarkerSize(0.7);

  h[0]->SetMarkerColor(kRed);
  h[1]->SetMarkerColor(kBlue);

  h[0]->SetMarkerStyle(22);
  h[1]->SetMarkerStyle(22);

  h[1]->Draw("hist");
  h[0]->Draw("hist sames");
  */

  hdQds[0]->Write();
  hdQds[1]->Write();

  hRatio[0]->Write();
  hRatio[1]->Write();

  //hRatio2D[0]->Write();
  //hRatio2D[1]->Write();

  dummy_file->Close();

  return;
}

//////////////////////////////// MAIN //////////////////////////////////////////

void dodqdx(string dirname = Dodqdx_Input, int run = 840, bool make_file = true ){
  //return the overall fitted gain distribution in root file for both views

  string filename = dirname + "select_tracks_"+to_string(run)+"_new_mitigation.root";
  
  bool optFit = false;

  cout << "optFit: " << optFit << endl;

  gStyle->SetOptFit(11111);
  gStyle->SetPalette(kRainBow);

  hdQds[0] = new TH1D("dQds_view_0", ";dQ/ds (fC/cm)", 200, 0, 100);
  hdQds[1] = new TH1D("dQds_view_1", ";dQ/ds (fC/cm)", 200, 0, 100);

  hRatio[0] = new TH1D("Ratio_view_0", ";(Hit Peak Amplitude)/(Hit Integral)", 200, -1, 1);
  hRatio[1] = new TH1D("Ratio_view_1", ";(Hit Peak Amplitude)/(Hit Integral)", 200, 1, 1);

  hRatio2D[0] = new TH2D("2D_Ratio_view_0", ";Hit Peak Amplitude; Hit Integral", 200, 0, 200, 200, 0, 2500);
  hRatio2D[1] = new TH2D("2D_Ratio_view_1", ";Hit Peak Amplitude; Hit Integral", 200, 0, 200, 200, 0, 2500);

  hRatioTheta[0]= new TProfile("RatioTheta_view_0", "View 0;#theta (rad);(Hit Peak Amplitude)/(Hit Integral)", 50, 0, 3.1415926/2);
  hRatioTheta[1]= new TProfile("RatioTheta_view_1", "View 1;#theta (rad);(Hit Peak Amplitude)/(Hit Integral)", 50, 0, 3.1415926/2);

  track *mip = new track();
  
  TFile *file = TFile::Open(filename.c_str());
  if(!file){ cout << "cannot open file: " << filename << endl; return; }

  TTree *tree = (TTree*)file->Get("mip_cut_delta");
  tree->SetBranchAddress("track", &mip);

  for(int ii=0; ii<tree->GetEntries(); ii++){
    tree->GetEntry(ii);

    for(int i =0; i< (int)mip->hits_trk.size(); i++){

      auto h = mip->hits_trk.at(i);
      //NB: free hits and hits aren't associated. Correspondence is 1:1
      //sorting migth be different
      auto free_h = mip->free_hits_trk.at(i);

      if ( isGood_lem(h.lem) ){

        if(h.dqdx >0)
          hdQds[h.view]->Fill(h.dqdx);

        if( isnormal(free_h.peak_amp/free_h.adc_integral) ){

          hRatio[free_h.view]->Fill(free_h.peak_amp/free_h.adc_integral);

          hRatio2D[free_h.view]->Fill(free_h.peak_amp, free_h.adc_integral);

          hRatioTheta[free_h.view]->Fill(mip->theta, free_h.peak_amp/free_h.adc_integral );
        }//adc integral protection
      }
    } //for hits
  }//for tree entries
  file->Close();

  // Setting fit range and start values
  double fr[2];
  double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  double chisqr;
  int    ndf;

  //////////////////////////////////////////////////////////////////////////////
  //quick fit with the overall dQds curve

  TF1 *function[NUM_OF_VIEWS];

  for(int view = 0; view < NUM_OF_VIEWS; view++ ){

    if(optFit){
      fr[0] = low[view];
      fr[1] = high[view];
    }
    else{
      FWHM(fr[0], fr[1], hdQds[view]);
    }

//    pllo[0]=0;
//    pllo[1]=0;
//    pllo[2]=1;
//    pllo[3]=0;
//    plhi[0]=10;
//    plhi[1]=20.0;
//    plhi[2]=1e+07;
//    plhi[3]=10;
//    sv[0]=1;
//    sv[1]=1;
//    sv[2]=1e+04;
//    sv[3]=1.0;

//    function_gain[view] = langaufit(hdQds[view],fr,sv,pllo,plhi,fp,fpe, &chisqr,&ndf, optFit);

    pllo[0]=0.1*hdQds[view]->GetStdDev();
    pllo[1]=0.1*hdQds[view]->GetMean();
    pllo[2]=0.1*hdQds[view]->Integral();
    pllo[3]=0.1*hdQds[view]->GetStdDev();
    plhi[0]=10*hdQds[view]->GetStdDev();
    plhi[1]=10*hdQds[view]->GetMean();
    plhi[2]=100*hdQds[view]->Integral();
    plhi[3]=10*hdQds[view]->GetStdDev();
    sv[0]=hdQds[view]->GetStdDev();
    sv[1]=hdQds[view]->GetMean();
    sv[2]=hdQds[view]->Integral();
    sv[3]=hdQds[view]->GetStdDev();
    function[view] = langaufit(hdQds[view],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf, false);

    cout << "view: " << view << " reduced chi^2: " << chisqr/ndf << endl;
  }

  if(make_file){
    string ofilename = Dodqdx_Output+"dqds_run_"+to_string(run)+"_new_mitigation.root";
    save_to_file(ofilename);
  }

}//end macro
