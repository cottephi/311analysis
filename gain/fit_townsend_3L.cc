#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>

//fit the 3L data with gain function. Parameters are transparency, A*rho and B*rho. LEM thickness is 1mm, variable is amp field in kV/cm
//Parameters found: T=1     A*rho=3907+/-11.45     B*rho=151.8+/-0.1107

void fit_townsend_3L(){
  
  vector<double> fields = {29,30,31,32,33,34,34.5,35};
  vector<double> fields_err = {};
  vector<double> gains = {7.8,11.7,18.5,30.4,51.4,89.6,122.1,166.3};
  vector<double> gains_errors = {};
  for(auto g : gains){
    gains_errors.push_back(g*0.02);
    fields_err.push_back(0);
  }
  TGraphErrors Gain_graph(fields.size(), &fields[0], &gains[0], &fields_err[0], &gains_errors[0]);
  TF1 *Theoretical_Gain = new TF1("Theoretical_Gain","[0]*exp([1]*0.1*exp(-[2]/x))",0,40);
  Theoretical_Gain->SetParameter(0,1);
  Theoretical_Gain->SetParName(0,"T");
  Theoretical_Gain->SetParLimits(0,0,1);
  Theoretical_Gain->SetParameter(1,4000);
  Theoretical_Gain->SetParName(1,"A*rho");
  Theoretical_Gain->SetParameter(2,200);
  Theoretical_Gain->SetParName(2,"B*rho");
  Gain_graph.Fit(Theoretical_Gain);
  Gain_graph.Draw("A*");
  gPad->SaveAs("./3L.png");
  return;
}
