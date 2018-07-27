////////////////////////////////////////////////////////////////////////////////
//
// Macro preparing the fft ch map for the DPhaseRealisticNoiseModule service
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
#include <TComplex.h>
#include <TFFTRealComplex.h>
#include <TFFTComplexReal.h>
#include <TVirtualFFT.h>

using namespace std;
using namespace art;

//global root variables
//TProfile2D *hFrequencyProfile2D;
//TProfile *hPowerSpectrum[2];
TProfile2D *hFrequencyProfile2D;
TProfile *hPowerSpectrum[2];

////////////////////////////// Graphics and i/o////////////////////////////

void load_style(){
  //to be improved: load custom style for the macro
  gStyle->SetPalette(kRainBow);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(11111);
}

void read_config(string &detector, int &n_views, int &ch_0, int &ch_1,
                                       float &sampling_freq, int &time_samples){

    /* Geometrical parameters set in a configuration file and read from it */

    detector = "3x1x1";
    n_views = 2;            // number of views
    ch_0 = 320;             //channels on view 0
    ch_1 = 960;             //channes on view 1
    sampling_freq = 2.5;    //sampling feq in MHz
    time_samples = 1667;    //number of time samples

    return;
  }

/////////////////////////////////Class definitiion /////////////////////////////

class fftUtils{
  public:
    fftUtils(int time_samples, float sampling_freq, bool roundup = false);
    ~fftUtils();

    //getters
    int   GetTimeSamples();
    float GetSamplingFreq();
    int   GetFrequencySize();

    void DoTimeArray(vector<float> &time_vector);
    void FFT(vector<float> input, vector<TComplex> & output);
    void PowerSpectrum(vector<TComplex> &fftdata, vector<float> &freq,
                                          vector<float> &spectrum, bool skipdc);
  private:
    int fDetectorTimeSize;
    int fSize;
    float fSampligFreq; //in MHz
    int fFreqSize;
    size_t nhalf;
    float df;

    TFFTRealComplex *fFFT;
};

/////////////////////////////  Class func //////////////////////////////////////

fftUtils::fftUtils(int time_samples, float sampling_freq, bool roundup){

  fDetectorTimeSize = time_samples;
  if(roundup)
    fSize = (int)( pow(2, ceil(log(fDetectorTimeSize)/log(2))) );
  else
    fSize = (int)( pow(2, floor(log(fDetectorTimeSize)/log(2))) );
  fSampligFreq = sampling_freq;
  fFreqSize = fSize/2+1;
  nhalf = fSize/2;
  df = fSampligFreq/fSize;

  fFFT = new TFFTRealComplex(fSize, false);
  int dummy[1] = {0};

  fFFT->Init("  ", -1, dummy);
}

fftUtils::~fftUtils(){ delete fFFT; }

int fftUtils::GetTimeSamples(){
  return fSize;
}

float fftUtils::GetSamplingFreq(){
  return fSampligFreq;
}

int fftUtils::GetFrequencySize(){
  return fFreqSize;
}

void fftUtils::DoTimeArray(vector<float> &time_vector){
  /* Fill up the time array up to match the number of fft points (for 311:
  1024 if roundup is false, 2048 if roundup is true  */

  int time_points = (int) time_vector.size();

  if(time_points > fSize){ return; } //nothing else to do here.

  // x2 mirror waveform ////////////////////////////////////////////////////////
  //fetch the first n_window entry of the array and calculate mean in the window

  int n_window = 200; double sum=0;
  for(int i = time_points-n_window; i<time_points; i++){ sum+=time_vector[i]; }
  double shift = sum/n_window; //shift for the mean within the window
  shift *=2;

  time_vector.resize(fSize);

  //do x2 mirror (default)
  for(int t = time_points; t<2*time_points; t++)
    time_vector[t] = -time_vector[2*time_points - t - 1 ] + shift;

  return;
}

void fftUtils::FFT(vector<float> input, vector<TComplex> & output){

    //Does the fft for the input at channel ch
    double real      = 0.;
    double imaginary = 0.;

    for(size_t p = 0; p < (size_t)fSize ; ++p)
      {
        if( p < (size_t)fSize )
          fFFT->SetPoint(p, input[p]);
        else // zero pad
          fFFT->SetPoint(p, 0);
      }

    fFFT->Transform();

    if( (int)output.size() < fFreqSize) output.resize( fFreqSize );

    for(int i = 0; i < fFreqSize; ++i){
      fFFT->GetPointComplex(i, real, imaginary);
      output[i] = TComplex(real, imaginary);
    }

    return;
}

void fftUtils::PowerSpectrum(vector<TComplex> &fftdata, vector<float> &freq,
                                          vector<float> &spectrum, bool skipdc){

  spectrum.clear();
  freq.clear();

  double p0 =0;
  size_t istart = 0;
  if(skipdc) istart++; //skip the first tdc value (sometimes does crazy things)

  for( size_t i=istart;i<fftdata.size();i++)
    {
      double p = fftdata[i].Rho();
      double f = df*i;
    //  if( i!=0 && i!=nhalf ) p *= 2;
      if( i==0 ) p0 = p;

      freq.push_back(f);
      spectrum.push_back(p);
    }
    return;
}

/////////////////////////////   Main    ////////////////////////////////////////

void dofft(string ifilename = "/eos/user/a/ascarpel/noise/raw/729-0_cnr32_cnr16.root ", string ofilename = "/eos/user/a/ascarpel/noise/models/fft_729-0_cnr32_cnr16.root", string inputtag = "caldata", unsigned int evt = 0){

  /* Do the fft transform for every channels, and fills up 2d histograms
  producing the frequency channle map and the average spectrum, save everything
  to root file. 0 means process all the events */

  //****************************************************************************
  
  string currentdir = get_current_dir_name();
  if (currentdir.find("pcotte") != std::string::npos){
    ifilename = "/eos/user/p/pcotte/311data/noise/raw/729-0_cnr32_cnr16.root";
    ofilepath = "/eos/user/p/pcotte/311data/noise/models/fft_729-0_cnr32_cnr16.root";
  }
  
  load_style();

  //assuming all these parameters read from config file in a second moment
  string detector;
  int n_views;
  int ch_0;       //channels on view 0
  int ch_1;       //channes on view 1
  float sampling_freq;  //sampling feq in MHz
  int time_samples;    //number of time samples

  read_config(detector, n_views, ch_0, ch_1, sampling_freq, time_samples);

  cout << "====== Dump detector info =======" << endl;
  cout << "Detector: " << detector << endl;
  cout << "Num of views: " << n_views << endl;
  cout << "Channels on view 0: " << ch_0 << endl;
  cout << "Channels on view 1 "  << ch_1 << endl;
  cout << "Sampling frequency (MHz): " << sampling_freq << endl;
  cout << "Number of time samples: " << time_samples << endl << endl;

  //fft
  fftUtils *fftUtil = new fftUtils(time_samples, sampling_freq, true);

  int fft_points = fftUtil->GetTimeSamples();
  int freq_size = fftUtil->GetFrequencySize();
  float df = sampling_freq/fft_points;

  cout << "====== Dump fft info =======" << endl;
  cout << "fft num of time samples: " << fft_points << endl;
  cout << "Max frequency (MHz): " << sampling_freq/2 << endl;
  cout << "Frequency size: " << freq_size << endl;
  cout << "Frequency resolution (MHz): " << df << endl << endl;

  //****************************************************************************
  //init histograms

  hFrequencyProfile2D = new TProfile2D( "hFrequencyProfile2D" ,
  "Channel frequency map; channels; f [MHz]" , ch_0+ch_1, 0, ch_0+ch_1, freq_size, 0,
                                                                 freq_size*df );

  hPowerSpectrum[0]  = new TProfile( "PowerSpectrum_0" , "View 0; f [MHz]",
                                                   freq_size, 0, freq_size*df );

  hPowerSpectrum[1]  = new TProfile( "PowerSpectrum_1" , "View 1; f [MHz]",
                                                   freq_size, 0, freq_size*df );

  TH1F *hfreq[freq_size];
  for(int f=0; f< freq_size; f++){
    string histname = "freq_"+to_string((int)f*df);
    hfreq[f] = new TH1F(histname.c_str(), histname.c_str(), 100, 0, 200);
  }

  //****************************************************************************
  //Begin the event-loop

  vector<float> time_vector;
  vector<TComplex> fftdata;
  vector<float> freq, spectrum;

  vector<string> filenames(1, ifilename);
  InputTag calwire_tag{ inputtag };

  unsigned int eventnum=0;
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {

    if(evt != 0){
      if( ev.eventAuxiliary().event() > evt){ break; }
      if(ev.eventAuxiliary().event() < evt){ ++eventnum; continue; }
    }

    auto const& calwires = *ev.getValidHandle<vector<recob::Wire>>(calwire_tag);

    if(calwires.empty()){
      cout << "Dataproduct calwire not found!" << endl;
      continue;
    }

    ++eventnum;
    cout << "Processing event: " << eventnum << endl;

    //begin the channel loop
    for(auto calwire : calwires){

      //re-initialize the vectors
      time_vector.clear();
      fftdata.clear();
      freq.clear(); spectrum.clear();

      time_vector = calwire.Signal();

      fftUtil->DoTimeArray(time_vector);

      fftUtil->FFT(time_vector, fftdata);
      fftUtil->PowerSpectrum( fftdata, freq, spectrum, false);

      int ch = calwire.Channel();

      //loop over frequecies and fill histograms
      for(int f=0; f<(int)freq.size(); f++){

          //skip here unwanted frequecies or spectrum values
          if(spectrum.at(f) == 0 || !isnormal(spectrum.at(f))){
            continue;
          } //skip 0s can be in dead ch.

          hfreq[f]->Fill( spectrum.at(f) );

          if (calwire.View() == 2){

            hFrequencyProfile2D->Fill( ch, freq.at(f), spectrum.at(f) );
            hPowerSpectrum[0]->Fill( freq.at(f), spectrum.at(f) );

          }
          else if(calwire.View() == 3){

            hFrequencyProfile2D->Fill( ch , freq.at(f), spectrum.at(f) );
            hPowerSpectrum[1]->Fill( freq.at(f), spectrum.at(f) );

          }
          else{
            cout << "not valid view " << endl;
          }

      }//end f
    } //end calwire
    cout << "Read next event " << endl;
  }//end event loop

  TCanvas *c = new TCanvas("c", "", 800, 600);
  c->cd(); c->SetLogz();
  hFrequencyProfile2D->Draw("colz");


  TCanvas *c1 = new TCanvas("c1", "", 1000, 900);
  c1->Divide(1,2);

  c1->cd(1);
  gPad->SetLogy();
  hPowerSpectrum[0]->Draw("hist");

  c1->cd(2);
  gPad->SetLogy();
  hPowerSpectrum[1]->Draw("hist");

  cout << "Save waveform to file " << ofilename << endl;
  TFile *fout = TFile::Open(ofilename.c_str(), "RECREATE");
  hPowerSpectrum[0]->Write();
  hPowerSpectrum[1]->Write();

  for(int f=0; f< freq_size; f++){
    hfreq[f]->Write();
  }

  fout->Close();

  cout << "All done!" << endl;
}//end macro
