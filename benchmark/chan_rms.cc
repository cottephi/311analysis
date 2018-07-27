////////////////////////////////////////////////////////////////////////////////
//
// Macro doing the chan rms and mean
// uses gallery
//
////////////////////////////////////////////////////////////////////////////////

// general header files:
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>

//gallery and larsoft
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "lardataobj/RecoBase/Wire.h"

// needed for ROOT routines:
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>

using namespace art;
using namespace std;

TH1F *hrms;
TH1F *hmean;
TH2F *hevd;

///////////////////////////////////////////////////////////////////////////////

void makemepretty(TH1F *h){

  h->SetMarkerStyle(20);
  h->SetMarkerSize(.4);
  h->SetMarkerColor(kBlue);
  h->SetLineColor(kBlue);
  h->SetStats(0);

  h->GetXaxis()->SetTitleFont(132);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.82);

  h->GetYaxis()->SetTitleFont(132);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.73);

}

///////////////////////////////////////////////////////////////////////////////

/*
double get_rms(vector<float> array){

  double sumxx=0.;
  double sumx=0.;
  double rms=0;

  for(auto adc : array){
    sumxx+=pow(adc,2);
    sumx+=adc;
  }

  cout << sumxx << " " << sumx*sumx << endl;

  double tval = sumxx-sumx*sumx;
  if(tval>0) rms = sqrt(tval/array.size());

  return rms;
}
*/

void GetMeanAndRMS(vector<float> array, double mean, double rms){

  //use welford method..consisted with qScan

  mean = 0;
  rms  = 0;

  double A = 0;
  double Q = 0;

  size_t istart = 3; //exclude first 3 point, can give crazy results

  for(size_t i=istart;i<array.size();i++){
    //cout<<data[i]<<endl;
    double d  = (double)array[i];
    double Ak = A + (d - A)/(i+1);
    double Qk = Q + (d - Ak)*(d-A);
    A = Ak;
    Q = Qk;
  }

    mean = A;
    rms  = sqrt( Q/(array.size()-1) );

    return;
}

///////////////////////////////////////////////////////////////////////////////

void chan_rms(string ifilename, string ofilename, size_t evt=0){
  /* Perform chan rms for all the channels in the geometry given with the input\
  TH2F files. evt flag allows to perform it only for a speficic event, if 0
  perform on all */

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  int ch0; int ch1; int tbins;
  ch0 = 320;
  ch1 = 960;
  tbins = 1667;

  string name = ";channel;RMS (ADC)";
  hrms = new TH1F("hrms", name.c_str(), ch0+ch1, -0.5, ch0+ch1-0.5);

  name = ";channel;Pedestal Mean (ADC)";
  hmean = new TH1F("hmean", name.c_str(), ch0+ch1, -0.5, ch0+ch1-0.5);

  vector<string> filenames(1, ifilename);

  InputTag calwire_tag{ "caldata" };

  unsigned int eventnum=0;

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {

    if(evt != 0){
      if( ev.eventAuxiliary().event() > evt){ break; }
      if(ev.eventAuxiliary().event() < evt){ ++eventnum; continue; }
    }

    cout << " Processing event " << ev.eventAuxiliary().event() << endl;

    eventnum++;

    auto const& calwires = *ev.getValidHandle<vector<recob::Wire>>(calwire_tag);

    if(calwires.empty()){
      cout << "Dataproduct calwire not found!" << endl;
      continue;
    }

    vector<float> time_array;
    int ch =0;
    double rms;
    double mean;

    for(auto calwire : calwires){

      time_array = calwire.Signal();

      GetMeanAndRMS(time_array, mean, rms);

      ch = calwire.Channel();

      hrms->Fill(ch, rms);
      hmean->Fill(ch, mean); //after ped subtraction in this case
    }

  }//end event loop

  if(evt == 0)
    hrms->Scale(1./eventnum); //normalize if calculated over many events

  if(evt == 0)
    hmean->Scale(1./eventnum); //normalize if calculated over many events

  makemepretty(hrms);
  hrms->Draw("hist p");

  TFile *ofile = TFile::Open(ofilename.c_str(), "RECREATE");
  hrms->Write();
  //hmean->Write();
  ofile->Close();

} //end macro
