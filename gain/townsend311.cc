////////////////////////////////////////////////////////////////////////////////
//
// Parametric fit of the 311 using townsend
//
////////////////////////////////////////////////////////////////////////////////

// general header files:
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
#include <THStack.h>

using namespace std;

TGraph *gQscan;
TGraph *g3L;

TF1 *fitf;

//////////////////////////////// Utilities /////////////////////////////////////

void draw(TCanvas *c, TGraph *g1, string header){

    gPad->SetTickx(1); gPad->SetTicky(1);

    g1->SetMarkerStyle(4);
    g1->SetMarkerSize(0.7);
    g1->GetXaxis()->SetTitle("LEM Field (kV/cm)");
    g1->GetYaxis()->SetTitle("LEM Amplification");
    g1->Draw("AP");
    c->Update();

    TF1 *f1 = g1->GetFunction("townsend");

    f1->SetLineWidth(2);
    f1->SetLineColor(kRed);
    f1->SetLineStyle(2);
    f1->Draw("LSAME");
    c->Update();

    TLegend *l = new TLegend(0.158635, 0.354667,  0.799197, 0.856);

    header = "\t"+header;
    l->SetHeader(header.c_str());
    l->AddEntry(g1, "Data", "p");
    l->AddEntry(f1, "Fit: ", "l");

    string label = "a="+to_string(f1->GetParameter(5));
    l->AddEntry((TObject*)0, label.c_str(), "");

    label = "b="+to_string(f1->GetParameter(6));
    l->AddEntry((TObject*)0, label.c_str(), "");

    label = "c="+to_string(f1->GetParameter(7));
    l->AddEntry((TObject*)0, label.c_str(), "");

    label = "l="+to_string(f1->GetParameter(1))+"+/-"+to_string(f1->GetParError(1));
    l->AddEntry((TObject*)0, label.c_str(), "");

    label = "T="+to_string(f1->GetParameter(0))+"+/-"+to_string(f1->GetParError(0));
    l->AddEntry((TObject*)0, label.c_str(), "");

    label = "#chi^2="+to_string(f1->GetChisquare());
    l->AddEntry((TObject*)0, label.c_str(), "");

    label = "ndf="+to_string(f1->GetNDF());
    l->AddEntry((TObject*)0, label.c_str(), "");

    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    c->Update();
}

void draw(TCanvas *c, TGraphErrors *g1, TGraph *g2){

    gPad->SetTickx(1); gPad->SetTicky(1);

    g1->SetMarkerStyle(4);
    g1->SetMarkerSize(0.7);
    g1->GetXaxis()->SetTitle("LEM Field (kV/cm)");
    g1->GetYaxis()->SetTitle("LEM Amplification");
    g1->Draw("AP");
    c->Update();

    g1->GetFunction("townsend")->SetLineWidth(2);
    g1->GetFunction("townsend")->SetLineColor(kRed);
    g1->GetFunction("townsend")->SetLineStyle(2);
    //gQscan->GetFunction("townsend")->SetRange(14, 38);
    g1->GetFunction("townsend")->Draw("same");
    c->Update();

    g2->SetLineColor(kBlack);
    g2->SetMarkerColor(kBlack);
    g2->SetMarkerSize(1.0);
    g2->GetXaxis()->SetLimits(22,38);
    g2->GetYaxis()->SetLimits(1,300);
    g2->GetXaxis()->SetTitle("LEM Field (kV/cm)");
    g2->GetYaxis()->SetTitle("LEM Amplification");
    g2->Draw("P");
    c->Update();

    g2->GetFunction("townsend")->SetLineWidth(2);
    g2->GetFunction("townsend")->SetLineColor(kBlue+2);
    g2->GetFunction("townsend")->SetLineStyle(3);
    //g3L->GetFunction("townsend")->SetRange(14, 38);
    g2->GetFunction("townsend")->Draw("same");
    c->Update();

    TLegend *l = new TLegend(0.158002, 0.59366,  0.505607, 0.851585);

    l->AddEntry(g1, "311 - data", "p");
    l->AddEntry(g1->GetFunction("townsend"), "311 - Fit (1.0 bar, 87 K)", "l");
    l->AddEntry(g2, "3L - data", "p" );
    l->AddEntry(g2->GetFunction("townsend"), "311 - Fit (0.980 bar, 87 K)", "l");
    l->Draw();

    c->Update();
}

void prepare_data(TGraph *g, double p){
  //TGraph *dummy = new TGraph();

  for(int i=0; i<g->GetN(); i++){
    double x=0.; double y=0.;
    g->GetPoint(i, x, y);
    g->SetPoint(i, x/p ,log(y)/p);
  }

  return;
}

/////////////////////////////// Functions  /////////////////////////////////////

double gainfunc(double *x, double *par){
  double e=x[0]*1000.;

  //fit parameters
  double T=par[0];
  double l=par[1];

  //fixed parameters
  double p=par[2];
  double t=par[3];
  double k=par[4];
  double a=par[5];
  double b=par[6];
  double c=par[7];

  return T*exp( a*p/t*exp(-b*p/t/(k*(e+c)) )*l);
}

double gainfunc_log(double *x, double *par){
  double e=x[0]*1000.;

  //fit parameters
  double T=par[0];
  double l=par[1];

  //fixed parameters
  double p=par[2];
  double t=par[3];
  double k=par[4];
  double a=par[5];
  double b=par[6];
  double c=par[7];

  return T + a*(p/t)*l*exp(-(b*(p/t))/(k*(e+c)) );
}

/////////////////////////////// Fetch data /////////////////////////////////////
TGraph *getDataQscan(){
  //fetch data from qScan
  cout << "Fetch qScan Lem gain data " << endl;
  string filename = "~/cernbox/field_scan/qscan_data.root";
  string currentdir = get_current_dir_name();
  if (currentdir.find("pcotte") != std::string::npos){
    filename = "/eos/user/p/pcotte/311analysis/gain/qscan_data.root";
  }

  TFile *f = new TFile(filename.c_str());

  TMultiGraph *tmpgraph = (TMultiGraph*)f->Get("qscan_gain");

  const int n =tmpgraph->GetListOfGraphs()->GetSize();

  //fetch the TMultiGraph to a normal TGraph
  TIter next( tmpgraph->GetListOfGraphs() );
  TGraph *g;
  vector<double> vx, vy;

  while( (g = (TGraph*)next()) ){
    for(int i=0; i<g->GetN(); i++){
      double x; double y;
      g->GetPoint(i, x, y);
      vx.push_back(x); vy.push_back(y);
    }
  }

  TGraph *gDummy = new TGraph();
  for(int i=0; i<vy.size(); i++){ gDummy->SetPoint(i, vx.at(i), vy.at(i) ); }

  return gDummy;
}

TGraphErrors *getDataLarsoft(){
  //fetch data from qScan
  cout << "Fetch LarSoft Lem gain data " << endl;

  string filename = "~/cernbox/field_scan/larsoft_gain.root ";
  string currentdir = get_current_dir_name();
  if (currentdir.find("pcotte") != std::string::npos){
    filename = "/eos/user/p/pcotte/311analysis/gain/larsoft_gain.root";
  }

  TFile *f = new TFile(filename.c_str());

  TGraphErrors *tmpgraph = (TGraphErrors*)f->Get("lem_larsoft");

  return tmpgraph;
}

TGraphErrors *getData3L(){
  //get the data from the 3L
  cout << "Fetch 3L Lem gain data " << endl;

  string filename = "~/cernbox/field_scan/summary.root";
  string currentdir = get_current_dir_name();
  if (currentdir.find("pcotte") != std::string::npos){
    filename = "/eos/user/p/pcotte/311analysis/gain/summary.root";
  }

  TFile *find = new TFile(filename.c_str());

  TGraphErrors *tmpgraph = new TGraphErrors();
  if(find->IsOpen()){
    tmpgraph = (TGraphErrors*)find->Get("lem_corrected");
  }else
  {
    cout << "cannot open file: " << filename << endl;
  }

  find->Close();
  cout << "Done" << endl;

  return tmpgraph;
}

void prepare_data(TGraph *g){
  //TGraph *dummy = new TGraph();

  for(int i=0; i<g->GetN(); i++){
    double x=0.; double y=0.;
    g->GetPoint(i, x, y);
    g->SetPoint( i, x ,log(y) );
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void study_c(double a, double b){

  gStyle->SetPalette(kRainBow);

  TCanvas *c_canvas = new TCanvas("c_canvas", "", 500, 400);
  gPad->SetTickx(1); gPad->SetTicky(2);

  TLegend *ll = new TLegend(0.1, 0.7,  0.3, 1.0);
  ll->SetNColumns(2);

  TF1 *func[200];

  //vary c in step of 100 from 0 to 10000
  for(int c=100; c<11000; c+=1000){

    cout << "function for c: " << c << endl;

    int ii = c/1000;

    //define the function
    func[ii] = new TF1("townsend", "gainfunc", 0, 40, 8);

    // fix parameters
    func[ii]->FixParameter(2, 1.0);  //pressure in bar
    func[ii]->FixParameter(3, 87);   //temprerature in K
    func[ii]->FixParameter(4, 0.95); //k (from the anode paper)

    //reduced a, b from the anode paper; c has to be estimated
    func[ii]->FixParameter(5, a);
    func[ii]->FixParameter(6, b);

    //set starting fit parameters
    func[ii]->SetParameter(0, 1.); //Transparency
    func[ii]->SetParameter(1, 0.07); //lem thichness (cm) from the anode paper

    func[ii]->FixParameter(7, c);

    func[ii]->SetLineColor(ii+1);
    func[ii]->GetYaxis()->SetTitle("Lem Aplification");
    func[ii]->GetXaxis()->SetTitle("E (kV/cm)");

    if(ii==0)
      func[ii]->Draw("L");
    else
      func[ii]->Draw("LSAME");

    cout << "f(0) " << func[ii]->Eval(0) << endl;

    string label = "c="+to_string(c);
    ll->AddEntry(func[ii], label.c_str(), "l");

    c_canvas->Update();
  }

  g3L->SetMarkerSize(0.7);
  g3L->SetMarkerColor(kBlack);
  g3L->SetMarkerStyle(20);
  g3L->Draw("P");

  ll->AddEntry(g3L, "3L - data", "p");
  ll->Draw();

  return;
}

///////////////////////////// //////////////////////////////////////////////////

void townsend311(string function="gainfunc_log"){

  //***************************************************************************
  //initial definitions of constants from MAGBOLZT

  //anode paiper
  double a = 274920;
  double b = 1.18668E+07;

  //lem paiper
  //double a = 651523;
  //double b = 1.62459E+07;

  //shu
  //double a=317461;
  //double b=1.2872448E+07;

  double c=0.;
  //double c=9379.;
  //double c = 7379.;

  //***************************************************************************

  //first of all fetch data

  if( function=="gainfunc" ){

    gQscan = new TGraph( *getDataQscan() );
    g3L = new TGraph( *getData3L() );

  }else if(function=="gainfunc_log"){

    gQscan = new TGraph( *getDataQscan() );
    g3L = new TGraph( *getData3L() );

    prepare_data( gQscan );
    prepare_data( g3L );

  }else{
    cout << "invalid option" << endl;
    return;
  }

  //study the fit function parametrisation
  //study_c(a, b);


  //define the function
  fitf = new TF1("townsend", function.c_str(), 14, 36, 8);

  // fix parameters
  fitf->FixParameter(2, 1.0);  //pressure in bar
  fitf->FixParameter(3, 87);   //temprerature in K
  fitf->FixParameter(4, 0.95); //k (from the anode paper)

  //reduced a, b from the anode paper; c has to be estimated
  fitf->FixParameter(5, a);
  fitf->FixParameter(6, b);
  fitf->FixParameter(7, c);

  //set starting fit parameters
  fitf->SetParameter(0, 1.); //Transparency
  fitf->SetParameter(1, 0.07); //lem thichness (cm) from the anode paper

  //now we fit the values of the 3L using c ~4100
  cout << "===== fit 3L data ======" << endl;

  fitf->SetRange(28, 36);
  fitf->FixParameter(7, 0); //after the study
  g3L->Fit("townsend", "R0");

  TCanvas *c3L = new TCanvas("c3L", "3L", 500, 400);
  draw(c3L, g3L, "3L");

  //now we fit the values of the 3L using c ~4100
  cout << "===== fit qScan data ======" << endl;

  fitf->SetRange(14, 32);
  fitf->FixParameter(7, 0); //after the study
  gQscan->Fit("townsend", "R0");

  TCanvas *c311 = new TCanvas("c311", "311", 500, 400);
  draw(c311, gQscan, "311");

}
