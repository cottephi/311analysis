#ifndef ROOT_TMyFileHeader
#define ROOT_TMyFileHeader


#ifndef ROOT_TObject
#include <TObject.h>
#endif
#include <TH1D.h>
#include <iostream>
//#include <stdlib.h>

class TMyFileHeader : public TObject {

protected:
  Int_t                                                        frun;
  Double_t                                                     fdrift;
  Double_t                                                     fextraction;
  Double_t                                                     famplification;
  Double_t                                                     finduction;
  Double_t                                                     frho;
  Double_t                                                     frho_var;
  Double_t                                                     fextr_eff_simu;
  Double_t                                                     ferr_extr_eff_simu;
  Double_t                                                     fextr_eff_gushin;
  Double_t                                                     ferr_extr_eff_gushin;
  Double_t                                                     find_eff;
  Double_t                                                     ferr_ind_eff;
  Double_t                                                     fdrift_correction_factor;
  Double_t                                                     ferr_drift_correction_factor;
  Double_t                                                     fdensity_correction_factor;
  Double_t                                                     ferr_density_correction_factor;
  Double_t                                                     ftotalcorrectionfactor_simu;
  Double_t                                                     ferrtotalcorrectionfactor_simu;
  Double_t                                                     ftotalcorrectionfactor_gushin;
  Double_t                                                     ferrtotalcorrectionfactor_gushin;
  Int_t                                                        fstart;
  Int_t                                                        fend;
  const char*                                                  fname;
  Double_t                                                     ftheoreticalgain;
  Double_t                                                     ftheoreticalgain3L;
  string                                                       ftrigger;
  map<string, pair<Double_t, Double_t> >                       fmpv;
  map<string, pair<Double_t, Double_t> >                       ferrmpv;
  map<string, map<int, pair<Double_t, Double_t> > >            fmpvlems;
  map<string, map<int, pair<Double_t, Double_t> > >            ferrmpvlems;
  map<string, Double_t>                                        fgaineff;
  map<string, Double_t>                                        ferrgaineff;
  map<string, map<int, Double_t> >                             fgainefflems;
  map<string, map<int, Double_t> >                             ferrgainefflems;
  map<string, Double_t>                                        fgainlem_simu;
  map<string, Double_t>                                        ferrgainlem_simu;
  map<string, map<int, Double_t> >                             fgainlem_simulems;
  map<string, map<int, Double_t> >                             ferrgainlem_simulems;
  map<string, Double_t >                                       fgainlem_gushin;
  map<string, Double_t >                                       ferrgainlem_gushin;
  map<string, map<int, Double_t> >                             fgainlem_gushinlems;
  map<string, map<int, Double_t> >                             ferrgainlem_gushinlems;
  map<string, Double_t >                                       ftransparancy;
  map<string, Double_t >                                       ferrtransparancy;
  map<string, map<int, Double_t> >                             ftransparancylems;
  map<string, map<int, Double_t> >                             ferrtransparancylems;
  map<string, pair<TH1D*, TH1D*> >                             fh;
  map<string, map<int, pair<TH1D*, TH1D*> > >                  fhlems;
  map<string, bool>                                            fgainprocessed;
  map<string, map<int, bool> >                                 fgainprocessedlems;
  map<string, map<pair<int,int>, pair<Double_t,Double_t> > >   fyzmpv;
  map<string, map<pair<int,int>, pair<Double_t,Double_t> > >   fyzerrmpv;
  map<string, map<pair<int,int>, pair<TH1D*,TH1D*> > >         fyzh;
  map<string, map<pair<int,int>, Double_t> >                   fyzgaineff;
  map<string, map<pair<int,int>, Double_t> >                   fyzerrgaineff;
  map<string, map<pair<int,int>, Double_t> >                   fyzgainlem_simu;
  map<string, map<pair<int,int>, Double_t> >                   fyzerrgainlem_simu;
  map<string, map<pair<int,int>, Double_t> >                   fyzgainlem_gushin;
  map<string, map<pair<int,int>, Double_t> >                   fyzerrgainlem_gushin;
  map<string, map<pair<int,int>, Double_t> >                   fyztransparancy;
  map<string, map<pair<int,int>, Double_t> >                   fyzerrtransparancy;
  map<string, map<pair<int,int>, bool> >                       fyzgainprocessed;

public:

  TMyFileHeader();
  TMyFileHeader(Int_t run, 
                Double_t drift, Double_t extraction, Double_t amplification, Double_t induction, 
                Int_t start = 0, Int_t end = 0, const char* name = 0, 
                Double_t rho = 0, Double_t rhovar = 0, 
                Double_t extr_eff_simu = -1, Double_t err_extr_eff_simu = -1,  Double_t extr_eff_gushin = -1,  Double_t err_extr_eff_gushin = -1, Double_t ind_eff = -1, Double_t err_ind_eff = -1, 
                Double_t drift_correction_factor = -1, Double_t err_drift_correction_factor = -1, Double_t density_correction_factor = -1, Double_t err_density_correction_factor = -1, Double_t theoreticalgain = -1, Double_t theoreticalgain3L = -1, string trigger = "NA");
  TMyFileHeader(const TMyFileHeader &header);
  virtual ~TMyFileHeader();

  Int_t                             GetRun()                   const {return frun;}
  Double_t                          GetDrift()                 const {return fdrift;}
  Double_t                          GetExtraction()            const {return fextraction;}
  Double_t                          GetAmplification()         const {return famplification;}
  Double_t                          GetInduction()             const {return finduction;}
  Int_t                             GetStartTime()             const {return fstart;}
  Int_t                             GetEndTime()               const {return fend;}
  Double_t                          GetRho()                   const {return frho;}
  Double_t                          GetRhoVar()                const {return frho_var;}
  Double_t                          GetExtrEffSimu()           const {return fextr_eff_simu;}
  Double_t                          GetErrExtrEffSimu()        const {return ferr_extr_eff_simu;}
  Double_t                          GetExtrEffGushin()         const {return fextr_eff_gushin;}
  Double_t                          GetErrExtrEffGushin()      const {return ferr_extr_eff_gushin;}
  Double_t                          GetIndEff()                const {return find_eff;}
  Double_t                          GetErrIndEff()             const {return ferr_ind_eff;}
  Double_t                          GetDriftCorr()             const {return fdrift_correction_factor;}
  Double_t                          GetErrDriftCorr()          const {return ferr_drift_correction_factor;}
  Double_t                          GetDensityCorr()           const {return fdensity_correction_factor;}
  Double_t                          GetErrDensityCorr()        const {return ferr_density_correction_factor;}
  Double_t                          GetTotalCorrSimu()         const {return ftotalcorrectionfactor_simu;}
  Double_t                          GetErrTotalCorrSimu()      const {return ferrtotalcorrectionfactor_simu;}
  Double_t                          GetTotalCorrGushin()       const {return ftotalcorrectionfactor_gushin;}
  Double_t                          GetErrTotalCorrGushin()    const {return ferrtotalcorrectionfactor_gushin;}
  Double_t                          GetTheoreticalGain()       const {return ftheoreticalgain;}
  Double_t                          GetTheoreticalGain3L()     const {return ftheoreticalgain3L;}
  string                            GetTrigger()               const {return ftrigger;}
  
  const char*                                       GetName()               const {return fname;}
  
  //////////////////////
  //Global get methods//
  //////////////////////
  
  map<string, pair<Double_t, Double_t> >                     GetMPVs()                                     const {return fmpv;}
  map<string, pair<Double_t, Double_t> >                     GetErrMPVs()                                  const {return ferrmpv;}
  map<string, map<int, pair<Double_t, Double_t> > >          GetMPVsLEMs()                                 const {return fmpvlems;}
  map<string, map<int, pair<Double_t, Double_t> > >          GetErrMPVsLEMs()                              const {return ferrmpvlems;}
  map<string, Double_t >                                     GetGainsEff()                                 const {return fgaineff;}
  map<string, Double_t >                                     GetErrGainsEff()                              const {return ferrgaineff;}
  map<string, map<int, Double_t > >                          GetGainsEffLEMs()                             const {return fgainefflems;}
  map<string, map<int, Double_t > >                          GetErrGainsEffLEMs()                          const {return ferrgainefflems;}
  map<string, Double_t >                                     GetGainsLemSimu()                             const {return fgainlem_simu;}
  map<string, Double_t >                                     GetErrGainsLemSimu()                          const {return ferrgainlem_simu;}
  map<string, map<int, Double_t > >                          GetGainsLemSimuLEMs()                         const {return fgainlem_simulems;}
  map<string, map<int, Double_t > >                          GetErrGainsLemSimuLEMs()                      const {return ferrgainlem_simulems;}
  map<string, Double_t >                                     GetGainsLemGushin()                           const {return fgainlem_gushin;}
  map<string, Double_t >                                     GetErrGainsLemGushin()                        const {return ferrgainlem_gushin;}
  map<string, map<int, Double_t > >                          GetGainsLemGushinLEMs()                       const {return fgainlem_gushinlems;}
  map<string, map<int, Double_t > >                          GetErrGainsLemGushinLEMs()                    const {return ferrgainlem_gushinlems;}
  map<string, Double_t >                                     GetTransparancies()                           const {return ftransparancy;}
  map<string, Double_t >                                     GetErrTransparancies()                        const {return ferrtransparancy;}
  map<string, map<int, Double_t > >                          GetTransparanciesLEMs()                       const {return ftransparancylems;}
  map<string, map<int, Double_t > >                          GetErrTransparanciesLEMs()                    const {return ferrtransparancylems;}
  map<string, pair<TH1D*, TH1D*> >                           GetHs()                                       const {return fh;}
  map<string, map<int, pair<TH1D*, TH1D*> > >                GetHsLEMs()                                   const {return fhlems;}
  map<string, bool>                                          AreGainsProcessed()                           const {return fgainprocessed;}
  map<string, map<int, bool> >                               AreGainsProcessedLEMs()                       const {return fgainprocessedlems;}
  map<string, map<pair<int,int>, pair<Double_t,Double_t> > > GetYZMPVs()                                   const {return fyzmpv;}
  map<string, map<pair<int,int>, pair<Double_t,Double_t> > > GetYZErrMPVs()                                const {return fyzerrmpv;}
  map<string, map<pair<int,int>, pair<TH1D*,TH1D*> > >       GetYZHs()                                     const {return fyzh;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZGainsEff()                               const {return fyzgaineff;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZErrGainsEff()                            const {return fyzerrgaineff;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZGainsLemSimu()                           const {return fyzgainlem_simu;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZErrGainsLemSimu()                        const {return fyzerrgainlem_simu;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZGainsLemGushin()                         const {return fyzgainlem_gushin;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZErrGainsLemGushin()                      const {return fyzerrgainlem_gushin;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZTransparancies()                         const {return fyztransparancy;}
  map<string, map<pair<int,int>, Double_t> >                 GetYZErrTransparancies()                      const {return fyzerrtransparancy;}
  map<string, map<pair<int,int>, bool> >                     AreYZGainsProcessed()                         const {return fyzgainprocessed;}
  
  
  ///////////////
  //set methods//
  ///////////////
  
  void SetRun(Int_t run)                                         {frun = run;                                                                                                    return;}
  void SetDrift(Double_t drift)                                  {fdrift = drift;                                                                                                return;}
  void SetExtraction(Double_t extraction)                        {fextraction = extraction;                                                                                      return;}
  void SetAmplification(Double_t amplification)                  {famplification = amplification;                                                                                return;}
  void SetInduction(Double_t induction)                          {finduction = induction;                                                                                        return;}
  void SetStartTime(Int_t t)                                     {fstart = t;                                                                                                    return;}
  void SetEndTime(Int_t t)                                       {fend = t;                                                                                                      return;}
  void SetRho(Double_t rho)                                      {frho = rho;                                                                                                    return;}
  void SetRhoVar(Double_t rhovar)                                {frho_var = rhovar;                                                                                             return;}
  void SetExtrEffSimu(Double_t extr_eff_simu)                    {fextr_eff_simu = extr_eff_simu;                                                                                return;}
  void SetErrExtrEffSimu(Double_t err_extr_eff_simu)             {ferr_extr_eff_simu = err_extr_eff_simu;                                                                        return;}
  void SetExtrEffGushin(Double_t extr_eff_gushin)                {fextr_eff_gushin = extr_eff_gushin;                                                                            return;}
  void SetErrExtrEffGushin(Double_t err_extr_eff_gushin)         {ferr_extr_eff_gushin = err_extr_eff_gushin;                                                                    return;}
  void SetIndEff(Double_t ind_eff)                               {find_eff = ind_eff;                                                                                            return;}
  void SetErrIndEff(Double_t err_ind_eff)                        {ferr_ind_eff = err_ind_eff;                                                                                    return;}
  void SetDriftCorr(Double_t drift_correction_factor)            {fdrift_correction_factor = drift_correction_factor;                                                            return;}
  void SetErrDriftCorr(Double_t err_drift_correction_factor)     {ferr_drift_correction_factor = err_drift_correction_factor;                                                    return;}
  void SetDensityCorr(Double_t density_correction_factor)        {fdensity_correction_factor = density_correction_factor;                                                        return;}
  void SetErrDensityCorr(Double_t err_density_correction_factor) {ferr_density_correction_factor = err_density_correction_factor;                                                return;}
  void SetTotalCorr_simu()                                       {ftotalcorrectionfactor_simu = fextr_eff_simu*find_eff*fdrift_correction_factor/fdensity_correction_factor;     return;}
  void SetErrTotalCorr_simu()                                    {ferrtotalcorrectionfactor_simu = pow(ferr_extr_eff_simu,2)/pow(fextr_eff_simu,2) + pow(ferr_ind_eff,2)/pow(find_eff,2) + pow(ferr_drift_correction_factor,2)/pow(fdrift_correction_factor,2) + pow(ferr_density_correction_factor,2)/pow(fdensity_correction_factor,2); return;}
  void SetTotalCorr_gushin()                                     {ftotalcorrectionfactor_gushin = fextr_eff_gushin*find_eff*fdrift_correction_factor/fdensity_correction_factor; return;}
  void SetErrTotalCorr_gushin()                                  {ferrtotalcorrectionfactor_gushin = pow(ferr_extr_eff_gushin,2)/pow(fextr_eff_gushin,2) + pow(ferr_ind_eff,2)/pow(find_eff,2) + pow(ferr_drift_correction_factor,2)/pow(fdrift_correction_factor,2) + pow(ferr_density_correction_factor,2)/pow(fdensity_correction_factor,2); return;}
  void SetName(const char* name)                                 {fname = name;                                                                                                  return;}
  void SetTheoreticalGain(Double_t theoreticalgain)              {ftheoreticalgain = theoreticalgain;                                                                            return;}
  void SetTheoreticalGain3L(Double_t theoreticalgain3L)          {ftheoreticalgain3L = theoreticalgain3L;                                                                        return;}
  void SetTrigger(string trigger)                                {ftrigger = trigger;                                                                                            return;}
  
  void SetMPV0(string cut_type, Double_t mpv){
    if(fmpv.find(cut_type) == fmpv.end()){fmpv[cut_type] = make_pair(-1,-1);}
    fmpv[cut_type].first = mpv;
    return;
  }
  void SetMPV1(string cut_type, Double_t mpv){
    if(fmpv.find(cut_type) == fmpv.end()){fmpv[cut_type] = make_pair(-1,-1);}
    fmpv[cut_type].second = mpv;
    return;
  }
  void SetErrMPV0(string cut_type, Double_t errmpv){
    if(ferrmpv.find(cut_type) == ferrmpv.end()){ferrmpv[cut_type] = make_pair(-1,-1);}
    ferrmpv[cut_type].first = errmpv;
    return;
  }
  void SetErrMPV1(string cut_type, Double_t errmpv){
    if(ferrmpv.find(cut_type) == ferrmpv.end()){ferrmpv[cut_type] = make_pair(-1,-1);}
    ferrmpv[cut_type].second = errmpv;
    return;
  }
  void SetMPVLEMs(string cut_type, map<int,pair<Double_t,Double_t> > mpvlems){
    fmpvlems[cut_type] = mpvlems;
    return;
  }
  void SetErrMPVLEMs(string cut_type, map<int,pair<Double_t,Double_t> > errmpvlems){
    ferrmpvlems[cut_type] = errmpvlems;
    return;
  }
  void SetMPVLEM(string cut_type, int lem, pair<Double_t,Double_t> mpvlem){
    fmpvlems[cut_type][lem] = mpvlem;
    return;
  }
  void SetErrMPVLEM(string cut_type, int lem, pair<Double_t,Double_t> errmpvlem){
    ferrmpvlems[cut_type][lem] = errmpvlem;
    return;
  }
  void SetH0(string cut_type, TH1D *h){
  
    if(fh.find(cut_type) == fh.end()){fh[cut_type] = make_pair((TH1D*)0,(TH1D*)0);}
    fh[cut_type].first = h;
    return;
  }
  void SetH1(string cut_type, TH1D *h){
    if(fh.find(cut_type) == fh.end()){fh[cut_type] = make_pair((TH1D*)0,(TH1D*)0);}
    fh[cut_type].second = h;
    return;
  }
  void SetHLEMs(string cut_type, map<int,pair<TH1D*,TH1D*> > hlems){
    fhlems[cut_type] = hlems;
    return;
  }
  void SetHLEM(string cut_type, int lem, pair<TH1D*,TH1D*> hlem){
    fhlems[cut_type][lem] = hlem;
    return;
  }
  void SetGainProcessed(string cut_type, bool processed){
    fgainprocessed[cut_type] = processed;
    return;
  }
  void SetGainProcessedLEMs(string cut_type, map<int, bool> processed){
    fgainprocessedlems[cut_type] = processed;
    return;
  }
  void SetGainProcessedLEM(string cut_type, int lem, bool processed){
    fgainprocessedlems[cut_type][lem] = processed;
    return;
  }
  void SetYZMPVs(string cut_type, map<pair<int,int>, pair<Double_t,Double_t> > yzmpvs){
    fyzmpv[cut_type] = yzmpvs;
    return;
  }
  void SetYZErrMPVs(string cut_type, map<pair<int,int>, pair<Double_t,Double_t> > yzerrmpvs){
    fyzerrmpv[cut_type] = yzerrmpvs;
    return;
  }
  void SetYZMPV(string cut_type, pair<int,int> YZ, pair<Double_t,Double_t> yzmpv){
    fyzmpv[cut_type][YZ] = yzmpv;
    return;
  }
  void SetYZErrMPV(string cut_type, pair<int,int> YZ, pair<Double_t,Double_t> yzerrmpv){
    fyzerrmpv[cut_type][YZ] = yzerrmpv;
    return;
  }
  void SetYZHs(string cut_type, map<pair<int,int>, pair<TH1D*,TH1D*> > yzhs){
    fyzh[cut_type] = yzhs;
    return;
  }
  void SetYZH(string cut_type, pair<int,int> YZ, pair<TH1D*,TH1D*> yzh){
    fyzh[cut_type][YZ] = yzh;
    return;
  }
  void SetYZGainProcesseds(string cut_type, map<pair<int,int>, bool > processeds){
    fyzgainprocessed[cut_type] = processeds;
    return;
  }
  void SetYZGainProcessed(string cut_type, pair<int,int> YZ, bool processed ){
    fyzgainprocessed[cut_type][YZ] = processed;
    return;
  }
  void ComputeYZGain(string cut_type, Double_t mpv_cos, /*Double_t err_mpv_cos,*/ int Y = 10000, int Z = 10000);
  
  void ComputeGain(string cut_type, Double_t mpv_cos, /*Double_t err_mpv_cos,*/ int lem = 0);
  
  void DrawHistos(string cut_type);
  void DrawHistosLEMs(string cut_type, int lem = 0);
  
  virtual void Print(Option_t *option="") const;
  void Clean();

  ClassDef(TMyFileHeader,1)
};

#endif
