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
  Int_t                                             frun;
  Double_t                                          fdrift;
  Double_t                                          fextraction;
  Double_t                                          famplification;
  Double_t                                          finduction;
  Double_t                                          frho;
  Double_t                                          frho_var;
  Double_t                                          fextr_eff_simu;
  Double_t                                          fextr_eff_gushin;
  Double_t                                          find_eff;
  Double_t                                          fdrift_correction_factor;
  Double_t                                          fdensity_correction_factor;
  Double_t                                          ftotalcorrectionfactor_simu;
  Double_t                                          ftotalcorrectionfactor_gushin;
  Int_t                                             fstart;
  Int_t                                             fend;
  const char*                                       fname;
  Double_t                                          ftheoreticalgain;
  Double_t                                          ftheoreticalgain3L;
  string                                            ftrigger;
  map<string, pair<Double_t, Double_t> >            fmpv;
  map<string, map<int, pair<Double_t, Double_t> > > fmpvlems;
  map<string, Double_t>                             fgaineff;
  map<string, map<int, Double_t> >                  fgainefflems;
  map<string, Double_t>                             fgainlem_simu;
  map<string, map<int, Double_t> >                  fgainlem_simulems;
  map<string, Double_t >                            fgainlem_gushin;
  map<string, map<int, Double_t> >                  fgainlem_gushinlems;
  map<string, Double_t >                            ftransparancy;
  map<string, map<int, Double_t> >                  ftransparancylems;
  map<string, pair<TH1D*, TH1D*> >                  fh;
  map<string, map<int, pair<TH1D*, TH1D*> > >       fhlems;
  map<string, bool>                                 fgainprocessed;
  map<string, map<int, bool> >                                 fgainprocessedlems;

public:

  TMyFileHeader();
  TMyFileHeader(Int_t run, 
                Double_t drift, Double_t extraction, Double_t amplification, Double_t induction, 
                Int_t start = 0, Int_t end = 0, const char* name = 0, 
                Double_t rho = 0, Double_t rhovar = 0, 
                Double_t extr_eff_simu = -1, Double_t extr_eff_gushin = -1, Double_t ind_eff = -1, 
                Double_t drift_correction_factor = -1, Double_t density_correction_factor = -1, Double_t theoreticalgain = -1, Double_t theoreticalgain3L = -1, string trigger = "NA");
  TMyFileHeader(const TMyFileHeader &header);
  virtual ~TMyFileHeader();

  Int_t                             GetRun()                const {return frun;}
  Double_t                          GetDrift()              const {return fdrift;}
  Double_t                          GetExtraction()         const {return fextraction;}
  Double_t                          GetAmplification()      const {return famplification;}
  Double_t                          GetInduction()          const {return finduction;}
  Int_t                             GetStartTime()          const {return fstart;}
  Int_t                             GetEndTime()            const {return fend;}
  Double_t                          GetRho()                const {return frho;}
  Double_t                          GetRhoVar()             const {return frho_var;}
  Double_t                          GetExtrEffSimu()        const {return fextr_eff_simu;}
  Double_t                          GetExtrEffGushin()      const {return fextr_eff_gushin;}
  Double_t                          GetIndEff()             const {return find_eff;}
  Double_t                          GetDriftCorr()          const {return fdrift_correction_factor;}
  Double_t                          GetDensityCorr()        const {return fdensity_correction_factor;}
  Double_t                          GetTotalCorrSimu()      const {return ftotalcorrectionfactor_simu;}
  Double_t                          GetTotalCorrGushin()    const {return ftotalcorrectionfactor_gushin;}
  Double_t                          GetTheoreticalGain()    const {return ftheoreticalgain;}
  Double_t                          GetTheoreticalGain3L()  const {return ftheoreticalgain3L;}
  string                            GetTrigger()            const {return ftrigger;}
  
  //////////////////////
  //Global get methods//
  //////////////////////
  
  map<string, pair<Double_t, Double_t> >            GetMPVs()               const {return fmpv;}
  map<string, map<int, pair<Double_t, Double_t> > > GetMPVsLEMs()           const {return fmpvlems;}
  map<string, Double_t >                            GetGainsEff()           const {return fgaineff;}
  map<string, map<int, Double_t > >                 GetGainsEffLEMs()       const {return fgainefflems;}
  map<string, Double_t >                            GetGainsLemSimu()       const {return fgainlem_simu;}
  map<string, map<int, Double_t > >                 GetGainsLemSimuLEMs()   const {return fgainlem_simulems;}
  map<string, Double_t >                            GetGainsLemGushin()     const {return fgainlem_gushin;}
  map<string, map<int, Double_t > >                 GetGainsLemGushinLEMs() const {return fgainlem_gushinlems;}
  map<string, Double_t >                            GetTransparancies()     const {return ftransparancy;}
  map<string, map<int, Double_t > >                 GetTransparanciesLEMs() const {return ftransparancylems;}
  map<string, pair<TH1D*, TH1D*> >                  GetHs()                 const {return fh;}
  map<string, map<int, pair<TH1D*, TH1D*> > >       GetHsLEMs()             const {return fhlems;}
  map<string, bool>                                 AreGainsProcessed()     const {return fgainprocessed;}
  map<string, map<int, bool> >                      AreGainsProcessedLEMs() const {return fgainprocessedlems;}
  
  const char*                                       GetName()               const {return fname;}
  
  ////////////////////////
  //cut_type get methods//
  ////////////////////////
  
  Double_t                          GetMPV0(string cut_type)              {
    if(fmpv.find(cut_type) == fmpv.end()){cout <<"No cut named " << cut_type << " found." << endl; return -1;}
    return fmpv[cut_type].first;
  }
  
  Double_t                          GetMPV1(string cut_type)              {
    if(fmpv.find(cut_type) == fmpv.end()){cout <<"No cut named " << cut_type << " found." << endl; return -1;}
    return fmpv[cut_type].second;
  }
  
  map<int,pair<Double_t,Double_t> > GetMPVLEMs(string cut_type)           {
    map<int,pair<Double_t,Double_t> > dummy;
    if(fmpvlems.find(cut_type) == fmpvlems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return fmpvlems[cut_type];
  }
  
  Double_t                          GetGainEff(string cut_type)           {
    if(fgaineff.find(cut_type) == fgaineff.end()){cout <<"No cut named " << cut_type << " found." << endl; return -1;}
    return fgaineff[cut_type];
  }
  
  map<int,Double_t>                 GetGainEffLEMs(string cut_type)       {
    map<int,Double_t> dummy;
    if(fgainefflems.find(cut_type) == fgainefflems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return fgainefflems[cut_type];
  }
  
  Double_t                          GetGainLemSimu(string cut_type)       {
    if(fgainlem_simu.find(cut_type) == fgainlem_simu.end()){cout <<"No cut named " << cut_type << " found." << endl; return -1;}
    return fgainlem_simu[cut_type];
  }
  
  map<int,Double_t>                 GetGainLemSimuLEMs(string cut_type)   {
    map<int,Double_t> dummy;
    if(fgainlem_simulems.find(cut_type) == fgainlem_simulems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return fgainlem_simulems[cut_type];
  }
  
  Double_t                          GetGainLemGushin(string cut_type)     {
    if(fgainlem_gushin.find(cut_type) == fgainlem_gushin.end()){cout <<"No cut named " << cut_type << " found." << endl; return -1;}
    return fgainlem_gushin[cut_type];
  }
  
  map<int,Double_t>                 GetGainLemGushinLEMs(string cut_type) {
    map<int,Double_t> dummy;
    if(fgainlem_gushinlems.find(cut_type) == fgainlem_gushinlems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return fgainlem_gushinlems[cut_type];
  }
  
  Double_t                          GetTransparancy(string cut_type)      {
    if(ftransparancy.find(cut_type) == ftransparancy.end()){cout <<"No cut named " << cut_type << " found." << endl; return -1;}
    return ftransparancy[cut_type];
  }
  
  map<int,Double_t>                 GetTransparancyLEMs(string cut_type)  {
    map<int,Double_t> dummy;
    if(ftransparancylems.find(cut_type) == ftransparancylems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return ftransparancylems[cut_type];
  }
  
  TH1D*                             GetH0(string cut_type)                {
    if(fh.find(cut_type) == fh.end()){cout <<"No cut named " << cut_type << " found." << endl; return NULL;}
    return fh[cut_type].first;
  }
  
  TH1D*                             GetH1(string cut_type)                {
    if(fh.find(cut_type) == fh.end()){cout <<"No cut named " << cut_type << " found." << endl; return NULL;}
    return fh[cut_type].second;
  }
  
  map<int,pair<TH1D*,TH1D*> >       GetHLEMs(string cut_type)             {
    map<int,pair<TH1D*,TH1D*> > dummy;
    if(fhlems.find(cut_type) == fhlems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return fhlems[cut_type];
  }
  
  bool                              IsGainProcessed(string cut_type)      {
    if(fgainprocessed.find(cut_type) == fgainprocessed.end()){cout <<"No cut named " << cut_type << " found." << endl; return false;}
    return fgainprocessed[cut_type];
  }
  
  map<int, bool>                    IsGainProcessedLEMs(string cut_type)  {
    map<int, bool> dummy;
    if(fgainprocessedlems.find(cut_type) == fgainprocessedlems.end()){cout <<"No cut named " << cut_type << " found." << endl; return dummy;}
    return fgainprocessedlems[cut_type];
  }
  
  ///////////////
  //set methods//
  ///////////////
  
  void SetRun(Int_t run)                                     {frun = run;                                                                                                    return;}
  void SetDrift(Double_t drift)                              {fdrift = drift;                                                                                                return;}
  void SetExtraction(Double_t extraction)                    {fextraction = extraction;                                                                                      return;}
  void SetAmplification(Double_t amplification)              {famplification = amplification;                                                                                return;}
  void SetInduction(Double_t induction)                      {finduction = induction;                                                                                        return;}
  void SetStartTime(Int_t t)                                 {fstart = t;                                                                                                    return;}
  void SetEndTime(Int_t t)                                   {fend = t;                                                                                                      return;}
  void SetRho(Double_t rho)                                  {frho = rho;                                                                                                    return;}
  void SetRhoVar(Double_t rhovar)                            {frho_var = rhovar;                                                                                             return;}
  void SetExtrEffSimu(Double_t extr_eff_simu)                {fextr_eff_simu = extr_eff_simu;                                                                                return;}
  void SetExtrEffGushin(Double_t extr_eff_gushin)            {fextr_eff_gushin = extr_eff_gushin;                                                                            return;}
  void SetIndEff(Double_t ind_eff)                           {find_eff = ind_eff;                                                                                            return;}
  void SetDriftCorr(Double_t drift_correction_factor)        {fdrift_correction_factor = drift_correction_factor;                                                            return;}
  void SetDensityCorr(Double_t density_correction_factor)    {fdensity_correction_factor = density_correction_factor;                                                        return;}
  void SetTotalCorr_simu()                                   {ftotalcorrectionfactor_simu = fextr_eff_simu*find_eff*fdrift_correction_factor/fdensity_correction_factor;     return;}
  void SetTotalCorr_gushin()                                 {ftotalcorrectionfactor_gushin = fextr_eff_gushin*find_eff*fdrift_correction_factor/fdensity_correction_factor; return;}
  void SetName(const char* name)                             {fname = name;                                                                                                  return;}
  void SetTheoreticalGain(Double_t theoreticalgain)          {ftheoreticalgain = theoreticalgain;                                                                            return;}
  void SetTheoreticalGain3L(Double_t theoreticalgain3L)      {ftheoreticalgain3L = theoreticalgain3L;                                                                        return;}
  void SetTrigger(string trigger)                            {ftrigger = trigger;                                                                                            return;}
  
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
  void SetMPVLEMs(string cut_type, map<int,pair<Double_t,Double_t> > mpvlems){
    fmpvlems[cut_type] = mpvlems;
    return;
  }
  void SetMPVLEM(string cut_type, int lem, pair<Double_t,Double_t> mpvlem){
    fmpvlems[cut_type][lem] = mpvlem;
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
  
  void ComputeGain(string cut_type, Double_t mpv_cos, int lem = 0);
  
  void DrawHistos(string cut_type);
  void DrawHistosLEMs(string cut_type, int lem = 0);
  
  virtual void Print(Option_t *option="") const;

  ClassDef(TMyFileHeader,1)
};

#endif
