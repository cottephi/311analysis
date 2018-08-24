#include "TClass.h"
#include "/eos/user/p/pcotte/311analysis/lib/TMyFileHeader.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TStyle.h>

ClassImp(TMyFileHeader)

TMyFileHeader::TMyFileHeader(): TObject(){
  frun                             = -1;
  fdrift                           = -1;
  fextraction                      = -1;
  famplification                   = -1;
  finduction                       = -1;
  fstart                           = 0;
  fend                             = 0;
  fname                            = 0;
  frho                             = 0;
  frho_var                         = 0;
  fextr_eff_simu                   = -1;
  ferr_extr_eff_simu               = 0;
  fextr_eff_gushin                 = -1;
  ferr_extr_eff_gushin             = 0;
  find_eff                         = -1;
  ferr_ind_eff                     = 0;
  fdrift_correction_factor         = -1;
  ferr_drift_correction_factor     = 0;
  fdensity_correction_factor       = -1;
  ferr_density_correction_factor   = 0;
  ftotalcorrectionfactor_simu      = -1;
  ferrtotalcorrectionfactor_simu  = 0;
  ftotalcorrectionfactor_gushin    = -1;
  ferrtotalcorrectionfactor_gushin = 0;
  ftheoreticalgain                 = -1;
  ftheoreticalgain3L               = -1;
  ftrigger                         = "NA";
  return;
}

TMyFileHeader::TMyFileHeader(Int_t run, Double_t drift, Double_t extraction, Double_t amplification, Double_t induction, Int_t start, Int_t end, const char* name, Double_t rho, Double_t rhovar, Double_t extr_eff_simu, Double_t err_extr_eff_simu, Double_t extr_eff_gushin, Double_t err_extr_eff_gushin, Double_t ind_eff, Double_t err_ind_eff, Double_t drift_correction_factor, Double_t err_drift_correction_factor, Double_t density_correction_factor, Double_t err_density_correction_factor, Double_t theoreticalgain, Double_t theoreticalgain3L, string trigger): TObject(){
  frun                             = run;
  fdrift                           = drift;
  fextraction                      = extraction;
  famplification                   = amplification;
  finduction                       = induction;
  fstart                           = start;
  fend                             = end;
  frho                             = rho;
  frho_var                         = rhovar;
  fextr_eff_simu                   = extr_eff_simu;
  ferr_extr_eff_simu               = err_extr_eff_simu;
  fextr_eff_gushin                 = extr_eff_gushin;
  ferr_extr_eff_gushin             = err_extr_eff_gushin;
  find_eff                         = ind_eff;
  ferr_ind_eff                     = err_ind_eff;
  fdrift_correction_factor         = drift_correction_factor;
  ferr_drift_correction_factor     = err_drift_correction_factor;
  fdensity_correction_factor       = density_correction_factor;
  ferr_density_correction_factor   = err_density_correction_factor;
  if(name == 0){fname              = ("run_"+to_string(run)).data();}
  else{fname                       = name;}
  ftotalcorrectionfactor_simu      = extr_eff_simu*ind_eff*drift_correction_factor/density_correction_factor;
  SetErrTotalCorr_simu();
  ftotalcorrectionfactor_gushin    = extr_eff_gushin*ind_eff*drift_correction_factor/density_correction_factor;
  SetErrTotalCorr_gushin();
  ftheoreticalgain                 = theoreticalgain;
  ftheoreticalgain3L               = theoreticalgain3L;
  ftrigger                         = trigger;
  return;
}

TMyFileHeader::TMyFileHeader(const TMyFileHeader &header): TObject(){
  frun                             = header.GetRun();
  fdrift                           = header.GetDrift();
  fextraction                      = header.GetExtraction();
  famplification                   = header.GetAmplification();
  finduction                       = header.GetInduction();
  fstart                           = header.GetStartTime();
  fend                             = header.GetEndTime();
  frho                             = header.GetRho();
  frho_var                         = header.GetRhoVar();
  fextr_eff_simu                   = header.GetExtrEffSimu();
  ferr_extr_eff_simu               = header.GetErrExtrEffSimu();
  fextr_eff_gushin                 = header.GetExtrEffGushin();
  ferr_extr_eff_gushin             = header.GetErrExtrEffGushin();
  find_eff                         = header.GetIndEff();
  ferr_ind_eff                     = header.GetErrIndEff();
  fdrift_correction_factor         = header.GetDriftCorr();
  ferr_drift_correction_factor     = header.GetErrDriftCorr();
  fdensity_correction_factor       = header.GetDensityCorr();
  ferr_density_correction_factor   = header.GetErrDensityCorr();
  fname                            = header.GetName();
  ftotalcorrectionfactor_simu      = header.GetTotalCorrSimu();
  ferrtotalcorrectionfactor_simu   = header.GetErrTotalCorrSimu();
  ftotalcorrectionfactor_gushin    = header.GetTotalCorrGushin();
  ferrtotalcorrectionfactor_gushin = header.GetErrTotalCorrGushin();
  ftheoreticalgain                 = header.GetTheoreticalGain();
  ftheoreticalgain3L               = header.GetTheoreticalGain3L();
  fmpv                             = header.GetMPVs();
  ferrmpv                          = header.GetErrMPVs();
  fmpvlems                         = header.GetMPVsLEMs();
  ferrmpvlems                      = header.GetErrMPVsLEMs();
  fgaineff                         = header.GetGainsEff();
  fgainefflems                     = header.GetGainsEffLEMs();
  fgainlem_simu                    = header.GetGainsLemSimu();
  fgainlem_simulems                = header.GetGainsLemSimuLEMs();
  fgainlem_gushin                  = header.GetGainsLemGushin();
  fgainlem_gushinlems              = header.GetGainsLemGushinLEMs();
  ftransparancy                    = header.GetTransparancies();
  ftransparancylems                = header.GetTransparanciesLEMs();
  fh                               = header.GetHs();
  fhlems                           = header.GetHsLEMs();
  fgainprocessed                   = header.AreGainsProcessed();
  fgainprocessedlems               = header.AreGainsProcessedLEMs();
  ftrigger                         = header.GetTrigger();
  fyzmpv                           = header.GetYZMPVs();
  fyzerrmpv                        = header.GetYZErrMPVs();
  fyzh                             = header.GetYZHs();
  fyzgaineff                       = header.GetYZGainsEff();
  fyzgainlem_simu                  = header.GetYZGainsLemSimu();
  fyzgainlem_gushin                = header.GetYZGainsLemGushin();
  fyztransparancy                  = header.GetYZTransparancies();
  fyzgainprocessed                 = header.AreYZGainsProcessed();
  return;
}

TMyFileHeader::~TMyFileHeader(){}

void TMyFileHeader::Print(Option_t *) const{
  std::cout << "----------- Run " << frun << " -----------" << std::endl;
  std::cout << "----------- Trigger " << ftrigger << " -----------" << std::endl;
  if(fstart == 0 || fend == 0){
    std::cout << "time not specified" << std::endl;
  }
  else{
    std::cout << "time: " << fstart << " to " << fend << " (duration: " << fstart-fend << " seconds)" << std::endl;
  }
  std::cout << "drift field: " << fdrift << " drift correction: " << fdrift_correction_factor << std::endl;
  std::cout << "extraction field: " << fextraction << " simulated extraction efficiency: " << fextr_eff_simu << " Gushin's extraction efficiency: " << fextr_eff_gushin <<  std::endl;
  std::cout << "amplification field: " << famplification << std::endl;
  std::cout << "induction field: " << finduction << " simulated induction efficiency: " << find_eff << std::endl;
  std::cout << "Rho: " << frho << " density correction factor: " << fdensity_correction_factor << std::endl;
  std::cout << "Rho Var: " << frho_var << std::endl;
  std::cout << "Total correction factor using simu extr eff: " << ftotalcorrectionfactor_simu << std::endl;
  std::cout << "Total correction factor using Gushin extr eff: " << ftotalcorrectionfactor_gushin << std::endl;
  std::cout << "Theoretical gain in LEM: " << ftheoreticalgain << std::endl;
  std::cout << "Theoretical gain in LEM at 3L's density: " << ftheoreticalgain3L << std::endl;
  std::cout << "------------------------------ " << std::endl;
  return;
}

void TMyFileHeader::ComputeGain(string cut_type, Double_t mpv_cos, /*errmpv_cos, */int lem){
  if(fmpv.find(cut_type) == fmpv.end()){
    cout <<"ComputeGain: No cut named " << cut_type << " found." << endl;
    return;
  }
  if(lem == 0){
    fgaineff[cut_type]                   = (fmpv[cut_type].first+fmpv[cut_type].second)/mpv_cos;
//    ferrgaineff[cut_type]                = TMath::Sqrt( (pow(ferrmpv[cut_type].first,2)+pow(fmpv[cut_type].second,2))/(pow(fmpv[cut_type].first+fmpv[cut_type].second),2) + pow(errmpv_cos,2)/pow(mpv_cos,2) );
    fgainlem_simu[cut_type]              = fgaineff[cut_type]/ftotalcorrectionfactor_simu;
//    ferrgainlem_simu[cut_type]                = TMath::Sqrt(pow(ferrgaineff[cut_type],2)/pow(fgaineff[cut_type],2) + pow(ferrtotalcorrectionfactor_simu,2)/pow(ftotalcorrectionfactor_simu,2));
    fgainlem_gushin[cut_type]            = fgaineff[cut_type]/ftotalcorrectionfactor_gushin;
//    ferrgainlem_gushin[cut_type]                = TMath::Sqrt(pow(ferrgaineff[cut_type],2)/pow(fgaineff[cut_type],2) + pow(ferrtotalcorrectionfactor_gushin,2)/pow(ftotalcorrectionfactor_gushin,2));
    ftransparancy[cut_type]              = fgaineff[cut_type]/ftheoreticalgain;
//    ferrtransparancy[cut_type]                = ferrgaineff[cut_type]/ftheoreticalgain;
  }
  else{
    if(fmpvlems.find(cut_type) == fmpvlems.end()){
      cout <<"ComputeGain: No cut named " << cut_type << " found in fmpvlems." << endl;
      return;
    }
    if(fmpvlems[cut_type].find(lem) == fmpvlems[cut_type].end()){
      cout <<"ComputeGain: LEM " << lem << " not found in cut_type " << cut_type << " found." << endl;
      return;
    }
    cout << "Computing gains of lEM " << lem << " for header." << endl;
    double mpv0                        = fmpvlems[cut_type][lem].first;
    double mpv1                        = fmpvlems[cut_type][lem].second;
    fgainefflems[cut_type][lem]        = (mpv0+mpv1)/mpv_cos;
//    ferrgainefflems[cut_type][lem]                = TMath::Sqrt( (pow(ferrmpvlems[cut_type][lem].first,2)+pow(ferrmpvlems[cut_type][lem].second,2))/pow((mpv0+mpv1),2) + pow(errmpv_cos,2)/pow(mpv_cos,2) );
    fgainlem_simulems[cut_type][lem]   = fgainefflems[cut_type][lem]/ftotalcorrectionfactor_simu;
//    ferrgainlem_simulems[cut_type]                = TMath::Sqrt(pow(ferrgainefflems[cut_type][lem],2)/pow(fgainefflems[cut_type][lem],2) + pow(ferrtotalcorrectionfactor_simu,2)/pow(ftotalcorrectionfactor_simu,2));
    fgainlem_gushinlems[cut_type][lem] = fgainefflems[cut_type][lem]/ftotalcorrectionfactor_gushin;
//    ferrgainlem_gushinlems[cut_type]                = TMath::Sqrt(pow(ferrgainefflems[cut_type][lem],2)/pow(fgainefflems[cut_type][lem],2) + pow(ferrtotalcorrectionfactor_gushin,2)/pow(ftotalcorrectionfactor_gushin,2));
    ftransparancylems[cut_type][lem]   = fgainefflems[cut_type][lem]/ftheoreticalgain;
//    ferrtransparancylems[cut_type]                = ferrgainefflems[cut_type][lem]/ftheoreticalgain;
  }
  return;
}

void TMyFileHeader::ComputeYZGain(string cut_type, Double_t mpv_cos, /*errmpv_cos, */int Y, int Z){
  if(fyzmpv.find(cut_type) == fyzmpv.end()){
    cout <<"ComputeYZGain: No cut named " << cut_type << " found." << endl;
    return;
  }
  pair<int,int> YZ = make_pair(Y,Z);
  if(Y != 10000 and Z != 10000){
    if(fyzmpv[cut_type].find(YZ) == fyzmpv[cut_type].end()){
      cout <<"ComputeYZGain: coordinates (" << Y << ";" << Z << ") not found in cut type " << cut_type << "." << endl;
      return;
    }
    cout << "Computing gains of (" << Y << ";" << Z << ") for header." << endl;
    double mpv0                        = fyzmpv[cut_type][YZ].first;
    double mpv1                        = fyzmpv[cut_type][YZ].second;
    fyzgaineff[cut_type][YZ]        = (mpv0+mpv1)/mpv_cos;
//    fyzerrgaineff[cut_type][lem]                = TMath::Sqrt( (pow(fyzerrmpv[cut_type][YZ].first,2)+pow(fyzerrmpv[cut_type][YZ].second,2))/pow((mpv0+mpv1),2) + pow(errmpv_cos,2)/pow(mpv_cos,2) );
    fyzgainlem_simu[cut_type][YZ]   = fyzgaineff[cut_type][YZ]/ftotalcorrectionfactor_simu;
//    fyzerrgainlem_simu[cut_type]                = TMath::Sqrt(pow(fyzerrgaineff[cut_type][YZ],2)/pow(fyzgaineff[cut_type][YZ],2) + pow(ferrtotalcorrectionfactor_simu,2)/pow(ftotalcorrectionfactor_simu,2));
    fyzgainlem_gushin[cut_type][YZ] = fyzgaineff[cut_type][YZ]/ftotalcorrectionfactor_gushin;
//    fyzerrgainlem_gushin[cut_type]                = TMath::Sqrt(pow(fyzerrgaineff[cut_type][YZ],2)/pow(fyzgaineff[cut_type][YZ],2) + pow(ferrtotalcorrectionfactor_gushin,2)/pow(ftotalcorrectionfactor_gushin,2));
    fyztransparancy[cut_type][YZ]   = fyzgaineff[cut_type][YZ]/ftheoreticalgain;
//    fyzerrtransparancy[cut_type]                = fyzerrgaineff[cut_type][YZ]/ftheoreticalgain;
  }
  else{
    cout << "ComputeYZGain: please specify both coordinates" << endl;
    return;
  }
  return;
}

void TMyFileHeader::DrawHistos(string cut_type){
  if(fh.find(cut_type) == fh.end()){
    cout <<"DrawHistos: No cut named " << cut_type << " found." << endl;
    return;
  }
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fh[cut_type].first->SetLineColor(kBlue);
  fh[cut_type].second->SetLineColor(kRed);
  double mymax = fh[cut_type].second->GetMaximum();
  if(fh[cut_type].first->GetMaximum() > mymax){mymax = fh[cut_type].first->GetMaximum();}
  fh[cut_type].second->SetTitle(string("run " + to_string(frun) + ";dQds(fC/cm);hits").data());
  fh[cut_type].first->SetTitle("view 0;dQds(fC/cm);hits");
  fh[cut_type].second->SetMaximum(1.1*mymax);
  fh[cut_type].first->Draw();
  fh[cut_type].second->Draw("SAME");
  return;
}

void TMyFileHeader::DrawHistosLEMs(string cut_type, int lem){
  if(fhlems.find(cut_type) == fhlems.end()){
    cout <<"DrawHistosLEMs: No cut named " << cut_type << " found." << endl;
    return;
  }
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  if(lem == 0){
    for(auto h : fhlems[cut_type]){
      h.second.first->SetLineColor(kBlue);
      h.second.second->SetLineColor(kRed);
      double max = h.second.second->GetMaximum();
      if(h.second.first->GetMaximum() > max){max = h.second.first->GetMaximum();}
      h.second.second->SetTitle(string("run " + to_string(frun) + " LEM " + to_string(h.first) + ";dQds(fC/cm);hits").data());
      h.second.first->SetTitle("view 0;dQds(fC/cm);hits");
      h.second.second->SetMaximum(1.1*max);
      h.second.second->Draw();
      h.second.first->Draw("SAME");
    }
  }
  else{
    if(fhlems[cut_type].find(lem) == fhlems[cut_type].end()){
      cout <<"DrawHistosLEMs: LEM  " << lem << " not found in cut type " << cut_type << "." << endl;
      return;
    }
    fhlems[cut_type][lem].first->SetLineColor(kBlue);
    fhlems[cut_type][lem].second->SetLineColor(kRed);
    double max = fhlems[cut_type][lem].second->GetMaximum();
    if(fhlems[cut_type][lem].first->GetMaximum() > max){max = fhlems[cut_type][lem].first->GetMaximum();}
    fhlems[cut_type][lem].second->SetTitle(string("run " + to_string(frun) + " LEM " + to_string(lem) + ";dQds(fC/cm);hits").data());
    fhlems[cut_type][lem].first->SetTitle("view 0;dQds(fC/cm);hits");
    fhlems[cut_type][lem].second->SetMaximum(1.1*max);
    fhlems[cut_type][lem].second->Draw();
    fhlems[cut_type][lem].first->Draw("SAME");
  }
  return;
}

void TMyFileHeader::Clean(){
  fmpv = {};
  ferrmpv = {};
  fmpvlems = {};
  ferrmpvlems = {};
  fgaineff = {};
  ferrgaineff = {};
  fgainefflems = {};
  ferrgainefflems = {};
  fgainlem_simu = {};
  ferrgainlem_simu = {};
  fgainlem_simulems = {};
  ferrgainlem_simulems = {};
  fgainlem_gushin = {};
  ferrgainlem_gushin = {};
  fgainlem_gushinlems = {};
  ferrgainlem_gushinlems = {};
  ftransparancy = {};
  ferrtransparancy = {};
  ftransparancylems = {};
  ferrtransparancylems = {};
  fh = {};
  fhlems = {};
  fgainprocessed = {};
  fgainprocessedlems = {};
  fyzmpv = {};
  fyzerrmpv = {};
  fyzh = {};
  fyzgaineff = {};
  fyzerrgaineff = {};
  fyzgainlem_simu = {};
  fyzerrgainlem_simu = {};
  fyzgainlem_gushin = {};
  fyzerrgainlem_gushin = {};
  fyztransparancy = {};
  fyzerrtransparancy = {};
  fyzgainprocessed = {};
  return;
}
