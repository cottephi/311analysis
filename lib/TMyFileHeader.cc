#include "TClass.h"
#include "/eos/user/p/pcotte/311analysis/lib/TMyFileHeader.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TStyle.h>

ClassImp(TMyFileHeader)

TMyFileHeader::TMyFileHeader(): TObject(){
  frun                          = -1;
  fdrift                        = -1;
  fextraction                   = -1;
  famplification                = -1;
  finduction                    = -1;
  fstart                        = 0;
  fend                          = 0;
  fname                         = 0;
  frho                          = 0;
  frho_var                      = 0;
  fextr_eff_simu                = -1;
  fextr_eff_gushin              = -1;
  find_eff                      = -1;
  fdrift_correction_factor      = -1;
  fdensity_correction_factor    = -1;
  ftotalcorrectionfactor_simu   = -1;
  ftotalcorrectionfactor_gushin = -1;
  ftheoreticalgain              = -1;
  ftheoreticalgain3L            = -1;
  ftrigger                      = "NA";
  return;
}

TMyFileHeader::TMyFileHeader(Int_t run, Double_t drift, Double_t extraction, Double_t amplification, Double_t induction, Int_t start, Int_t end, const char* name, Double_t rho, Double_t rhovar, Double_t extr_eff_simu, Double_t extr_eff_gushin, Double_t ind_eff, Double_t drift_correction_factor, Double_t density_correction_factor, Double_t theoreticalgain, Double_t theoreticalgain3L, string trigger): TObject(){
  frun                          = run;
  fdrift                        = drift;
  fextraction                   = extraction;
  famplification                = amplification;
  finduction                    = induction;
  fstart                        = start;
  fend                          = end;
  frho                          = rho;
  frho_var                      = rhovar;
  fextr_eff_simu                = extr_eff_simu;
  fextr_eff_gushin              = extr_eff_gushin;
  find_eff                      = ind_eff;
  fdrift_correction_factor      = drift_correction_factor;
  fdensity_correction_factor    = density_correction_factor;
  if(name == 0){fname           = ("run_"+to_string(run)).data();}
  else{fname                    = name;}
  ftotalcorrectionfactor_simu   = extr_eff_simu*ind_eff*drift_correction_factor/density_correction_factor;
  ftotalcorrectionfactor_gushin = extr_eff_gushin*ind_eff*drift_correction_factor/density_correction_factor;
  ftheoreticalgain              = theoreticalgain;
  ftheoreticalgain3L            = theoreticalgain3L;
  ftrigger                      = trigger;
  return;
}

TMyFileHeader::TMyFileHeader(const TMyFileHeader &header): TObject(){
  frun                          = header.GetRun();
  fdrift                        = header.GetDrift();
  fextraction                   = header.GetExtraction();
  famplification                = header.GetAmplification();
  finduction                    = header.GetInduction();
  fstart                        = header.GetStartTime();
  fend                          = header.GetEndTime();
  frho                          = header.GetRho();
  frho_var                      = header.GetRhoVar();
  fextr_eff_simu                = header.GetExtrEffSimu();
  fextr_eff_gushin              = header.GetExtrEffGushin();
  find_eff                      = header.GetIndEff();
  fdrift_correction_factor      = header.GetDriftCorr();
  fdensity_correction_factor    = header.GetDensityCorr();
  fname                         = header.GetName();
  ftotalcorrectionfactor_simu   = header.GetTotalCorrSimu();
  ftotalcorrectionfactor_gushin = header.GetTotalCorrGushin();
  ftheoreticalgain              = header.GetTheoreticalGain();
  ftheoreticalgain3L            = header.GetTheoreticalGain3L();
  fmpv                          = header.GetMPVs();
  fmpvlems                      = header.GetMPVsLEMs();
  fgaineff                      = header.GetGainsEff();
  fgainefflems                  = header.GetGainsEffLEMs();
  fgainlem_simu                 = header.GetGainsLemSimu();
  fgainlem_simulems             = header.GetGainsLemSimuLEMs();
  fgainlem_gushin               = header.GetGainsLemGushin();
  fgainlem_gushinlems           = header.GetGainsLemGushinLEMs();
  ftransparancy                 = header.GetTransparancies();
  ftransparancylems             = header.GetTransparanciesLEMs();
  fh                            = header.GetHs();
  fhlems                        = header.GetHsLEMs();
  fgainprocessed                = header.AreGainsProcessed();
  fgainprocessedlems            = header.AreGainsProcessedLEMs();
  ftrigger                      = header.GetTrigger();
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

void TMyFileHeader::ComputeGain(string cut_type, Double_t mpv_cos, int lem){
  if(fmpv.find(cut_type) == fmpv.end()){
    cout <<"No cut named " << cut_type << " found." << endl;
    return;
  }
  if(lem == 0){
    fgaineff[cut_type]                   = (fmpv[cut_type].first+fmpv[cut_type].second)/mpv_cos;
    fgainlem_simu[cut_type]              = fgaineff[cut_type]/ftotalcorrectionfactor_simu;
    fgainlem_gushin[cut_type]            = fgaineff[cut_type]/ftotalcorrectionfactor_gushin;
    ftransparancy[cut_type]              = fgaineff[cut_type]/ftheoreticalgain;
  }
  else{
    double mpv0                        = fmpvlems[cut_type][lem].first;
    double mpv1                        = fmpvlems[cut_type][lem].second;
    fgainefflems[cut_type][lem]        = (mpv0+mpv1)/mpv_cos;
    fgainlem_simulems[cut_type][lem]   = fgainefflems[cut_type][lem]/ftotalcorrectionfactor_simu;
    fgainlem_gushinlems[cut_type][lem] = fgainefflems[cut_type][lem]/ftotalcorrectionfactor_gushin;
    ftransparancylems[cut_type][lem]   = fgainefflems[cut_type][lem]/ftheoreticalgain;
  }
  return;
}

void TMyFileHeader::DrawHistos(string cut_type){
  if(fh.find(cut_type) == fh.end()){
    cout <<"No cut named " << cut_type << " found." << endl;
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
    cout <<"No cut named " << cut_type << " found." << endl;
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
      cout <<"LEM  " << lem << " not found." << endl;
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
