#ifndef ROOT_TMyFileHeader
#define ROOT_TMyFileHeader


#ifndef ROOT_TObject
#include <TObject.h>
#endif
#include <iostream>
//#include <stdlib.h>

class TMyFileHeader : public TObject {

protected:
  Int_t       frun;
  Double_t    fdrift;
  Double_t    fextraction;
  Double_t    famplification;
  Double_t    finduction;
  Double_t    frho;
  Double_t    frho_var;
  Int_t       fstart;
  Int_t       fend;
  const char* fname;

public:

  TMyFileHeader();
  TMyFileHeader(Int_t run, Double_t drift, Double_t extraction, Double_t amplification, Double_t induction, Int_t start = 0, Int_t end = 0, const char* name = 0, Double_t rho = 0, Double_t rhovar = 0);
  TMyFileHeader(const TMyFileHeader &header);
  virtual ~TMyFileHeader();

  Int_t                GetRun()           const {return frun;}
  Double_t             GetDrift()         const {return fdrift;}
  Double_t             GetExtraction()    const {return fextraction;}
  Double_t             GetAmplification() const {return famplification;}
  Double_t             GetInduction()     const {return finduction;}
  Int_t                GetStartTime()     const {return fstart;}
  Int_t                GetEndTime()       const {return fend;}
  Double_t             GetRho()           const {return frho;}
  Double_t             GetRhoVar()        const {return frho_var;}
  const char*          GetName()          const {return fname;}
  
  void                 SetRun(Int_t run)                        {frun = run;                     return;}
  void                 SetDrift(Double_t drift)                 {fdrift = drift;                 return;}
  void                 SetExtraction(Double_t extraction)       {fextraction = extraction;       return;}
  void                 SetAmplification(Double_t amplification) {famplification = amplification; return;}
  void                 SetInduction(Double_t induction)         {finduction = induction;         return;}
  void                 SetStartTime(Int_t t)                    {fstart = t;                     return;}
  void                 SetEndTime(Int_t t)                      {fend = t;                       return;}
  void                 SetRho(Double_t rho)                     {frho = rho;                     return;}
  void                 SetRhoVar(Double_t rhovar)               {frho_var = rhovar;              return;}
  void                 SetName(const char* name)                {fname = name;                   return;}

  virtual void         Print(Option_t *option="") const;

  ClassDef(TMyFileHeader,1)
};

#endif
