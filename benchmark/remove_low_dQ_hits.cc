////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}


void remove_low_dQ_hits(int run = 840, string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_dQ = m_dQ;
  method_ds = m_ds;
  if(!Load_Version(version)){return;}

  gErrorIgnoreLevel = kError;
//  gStyle->SetOptStat(0);
//  gStyle->SetStatY(0.35);                
//  // Set y-position (fraction of pad size)
//  gStyle->SetStatX(0.85);                
//  // Set x-position (fraction of pad size)
//  gStyle->SetStatW(0.2);                
//  // Set width of stat-box (fraction of pad size)
//  gStyle->SetStatH(0.15);                
  // Set height of stat-box (fraction of pad size)

  
  TH1D hdQds;
  string inputpath = dQds_Output + "before_cuts/" + to_string(run) + ".root";
  #if verbose
  cout << "Opening file " << inputpath << endl;
  #endif
  TFile ifile(inputpath.data(), "UPDATE");
  if(!ifile.IsOpen()){
    cout << "ERROR: can not open " << inputpath << endl;
    return;
  }
  TH1D *dQds = 0;
  
  #if verbose
  cout << "Removing previous corrected plots... " << endl;
  #endif
  
  for(auto obj : *(TList*)ifile.GetListOfKeys()){
    string name = ((TKey*)obj)->GetName();
    int i=2;
    cout << name << endl;
    TDirectory *dir = ifile.GetDirectory("/");
    while(i<10){
      string to_remove = name+";"+to_string(i);
      #if verbose
      cout << "Removing " << to_remove << endl;
      #endif
      dir->Delete(to_remove.data());
      i++;
    }
    i=0;
    dQds = 0; delete dQds;
  }

  for(auto obj : *(TList*)ifile.GetListOfKeys()){
    string name = ((TKey*)obj)->GetName();
    if(name.find("lowdqhitsremoved") != string::npos){
      int i=1;
      TDirectory *dir = ifile.GetDirectory("/");
      while(true){
        string to_remove = name+";"+to_string(i);
        #if verbose
        cout << "Removing " << to_remove << endl;
        #endif
        dir->Delete(to_remove.data());
        TH1D *dummy = 0;
        ifile.GetObject(name.data(),dummy);
        if(dummy == 0){break;}
        dummy = 0; delete dummy;
        i++;
      }
    }
    dQds = 0; delete dQds;
  }
  
  #if verbose
  cout << "Applying correction... " << endl;
  #endif
  for(auto obj : *(TList*)ifile.GetListOfKeys()){
    dQds = (TH1D*)((TKey*)obj)->ReadObj();
    if(string(dQds->GetName()).find("lowdqhitsremoved") != string::npos){
      dQds = 0; delete dQds;
      break;
    }
    hdQds = TH1D(*dQds);
    hdQds.SetDirectory(0);
    hdQds.SetName(string(string(hdQds.GetName())+"_lowdqhitsremoved").data());
    TF1 *myfit = new TF1("myfit","[0]*exp(-[1]*x)",0,20);
    myfit->SetParameter(0,12.);
    myfit->SetParameter(1,0.5);
    myfit->SetParLimits(1,-1000,0);
    dQds->Fit(myfit,"","",0,2);
    double pvalue = myfit->GetProb();
    if(pvalue < 0.05 or pvalue != pvalue){continue;}
    for(int i = 1; i<=hdQds.GetNbinsX(); i++){
      double x = hdQds.GetBinCenter(i);
      int bincontent = hdQds.GetBinContent(i);
      double y;
      double xmin, xmax;
      myfit->GetRange(xmin,xmax);
      if(x<xmin){continue;}
      if(x>xmax){break;}
      y = myfit->Eval(x);
      if(bincontent-y < 0){hdQds.SetBinContent(i,0);}
      else{hdQds.SetBinContent(i,bincontent-y);}
    }
    dQds->Write(dQds->GetName(),TObject::kOverwrite);
    hdQds.Write();
    dQds = 0; delete dQds;
  }
  ifile.Close();
  #if verbose 
  cout << "...done" << endl;
  #endif
  return;
}//end macro
