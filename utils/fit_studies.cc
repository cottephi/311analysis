#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <cmath>


#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <TStyle.h>

#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveText.h>
using namespace std;


#include <sstream>
using std::ostringstream;
#include <fstream>



Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1* langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[4]  reasonable start values for the fit
  //   parlimitslo[4]  lower parameter limits
  //   parlimitshi[4]  upper parameter limits
  //   fitparams[4]    returns the final fit parameters
  //   fiterrors[4]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf

  Int_t i;
  TString FunName = TString::Format("Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName.Data());
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName.Data(),langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");

  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  his->Fit(FunName.Data(),"RB0Q");   // fit within specified range, use ParLimits, do not plot

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf

  return (ffit);              // return fit function

}






void fit_studies(){

  double mpv= 11;
  double Lwidth= 0.95;
  double Gwidth = 1.6;

  double norm = 1;

  cout << " mpv = " << mpv << " Lwidth = " << Lwidth << " = " << Gwidth << endl;

  TH1D* h_mpv = new TH1D("h_mpv", "MPV; Pull", 200, -3, 3);
  TH1D* h_lwidth = new TH1D("h_lwidth", "L. Width; Pull", 200, -3, 3);
  TH1D* h_gwidth = new TH1D("h_gwidth", "G. Width; Pull", 200, -3, 3);


  TH1F* hin = new TH1F("hin", "; charge; ", 200, 0., 100);

  TF1* fin = new TF1("fin", langaufun, 0., 100., 4);
  fin->SetParameters(Lwidth, mpv, norm, Gwidth);




  Double_t range[2];
  Double_t start[4];
  Double_t minv[4];
  Double_t maxv[4];
  Double_t fitp[4];
  Double_t fite[4];
  Double_t chi2;
  int ndf;

  range[0] = 0.;
  range[1] = 25.;
  start[0] = 1.; start[1] = 2.; start[2] = 100.; start[3] = 1.;
  minv[0] = 0.; minv[1] = 0.5; minv[2] = 0.; minv[3] = 0.;
  maxv[0] = 10.; maxv[1] = 25.; maxv[2] = 20000.; maxv[3] = 10.;
  TF1* fout;

  for(int kiter = 0; kiter < 5000; kiter++){

    hin->FillRandom("fin", 10000);
    fout = langaufit(hin, range, start, minv, maxv, fitp, fite, &chi2, &ndf);
    h_mpv->Fill( (fout->GetParameter(1) - mpv)/fout->GetParError(1));
    h_lwidth->Fill( (fout->GetParameter(0) - Lwidth)/fout->GetParError(0));
    h_gwidth->Fill( (fout->GetParameter(3) - Gwidth)/fout->GetParError(3));
    hin->Reset();

    if(kiter%200==0 and kiter!=0) cout << kiter << " iteration done ! " << endl;
  }
  cout << " 5000 iteration done ! " << endl;

  TFile *fres = new TFile(Form("pulls_mpv_%.1f_lw_%.1f_gw_%.1f.root",mpv,Lwidth,Gwidth),"recreate");
  fres->cd();
  hin->FillRandom("fin", 10000);
  hin->Write();
  h_mpv->Write();
  h_lwidth->Write();
  h_gwidth->Write();
  fres->Close();
  return;

  /*
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);


  TCanvas *c = new TCanvas("c","c",1200, 450);
  //  hin->Draw();
  //  fout->Draw("same");
  c->Divide(3,1);
  c->cd(1);
  h_mpv->Draw();
  h_mpv->Fit("gaus");
  c->cd(2);
  h_lwidth->Draw();
  h_lwidth->Fit("gaus");
  c->cd(3);
  h_gwidth->Draw();
  h_gwidth->Fit("gaus");
  c->SaveAs("test.pdf");
*/
}
