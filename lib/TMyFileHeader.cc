#include "TClass.h"
#include "/eos/user/p/pcotte/311analysis/lib/TMyFileHeader.h"

ClassImp(TMyFileHeader)

TMyFileHeader::TMyFileHeader(): TObject(){
  frun = -1;
  fdrift = -1;
  fextraction = -1;
  famplification = -1;
  finduction = -1;
  fstart = 0;
  fend = 0;
  fname = 0;
  frho = 0;
  frho_var = 0;
  return;
}

TMyFileHeader::TMyFileHeader(Int_t run, Double_t drift, Double_t extraction, Double_t amplification, Double_t induction, Int_t start, Int_t end, const char* name, Double_t rho, Double_t rhovar): TObject(){
  frun = run;
  fdrift = drift;
  fextraction = extraction;
  famplification = amplification;
  finduction = induction;
  fstart = start;
  fend = end;
  frho = rho;
  frho_var = rhovar;
  if(name == 0){
    fname = ("run_"+to_string(run)).data();
  }
  else{
    fname = name;
  }
  return;
}

TMyFileHeader::TMyFileHeader(const TMyFileHeader &header): TObject(){
  frun = header.GetRun();
  fdrift = header.GetDrift();
  fextraction = header.GetExtraction();
  famplification = header.GetAmplification();
  finduction = header.GetInduction();
  fstart = header.GetStartTime();
  fend = header.GetEndTime();
  frho = header.GetRho();
  frho_var = header.GetRhoVar();
  fname = header.GetName();
  return;
}

TMyFileHeader::~TMyFileHeader(){}

void TMyFileHeader::Print(Option_t *) const{
  std::cout << "----------- Run " << frun << " -----------" << std::endl;
  if(fstart == 0 || fend == 0){
    std::cout << "time not specified" << std::endl;
  }
  else{
    std::cout << "time: " << fstart << " to " << fend << " (duration: " << fstart-fend << " seconds)" << std::endl;
  }
  std::cout << "drift field: " << fdrift << std::endl;
  std::cout << "extraction field: " << fextraction << std::endl;
  std::cout << "amplification field: " << famplification << std::endl;
  std::cout << "induction field: " << finduction << std::endl;
  std::cout << "Rho: " << frho << std::endl;
  std::cout << "Rho Var: " << frho_var << std::endl;
  std::cout << "------------------------------ " << std::endl;
  return;
}
