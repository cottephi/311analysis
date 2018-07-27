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
#include <TMyFileHeader.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}


void charging_up_all(vector<int> run_list = {771,840,842,993,996,998,1009,1011,1016,1035,1036,1037}, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  gErrorIgnoreLevel = kError;

  pair<vector<double>,vector<double> > slope_0;
  pair<vector<double>,vector<double> > intercept_0;
  pair<vector<double>,vector<double> > slope_1;
  pair<vector<double>,vector<double> > intercept_1;
  map<int, pair<vector<double>,vector<double> > > slope_0_ByLEMs;
  map<int, pair<vector<double>,vector<double> > > intercept_0_ByLEMs;
  map<int, pair<vector<double>,vector<double> > > slope_1_ByLEMs;
  map<int, pair<vector<double>,vector<double> > > intercept_1_ByLEMs;
  
  for(auto run : run_list){
  
    string header_file = runs_headers + to_string(run)+".root";
    if(!ExistTest(header_file)){gROOT->ProcessLine(string(".x " + path_311analysis + "utils/save_runs_headers.cc({" + to_string(run) + "})").data());}
    TFile ifile_header(header_file.data(),"READ");
    TMyFileHeader *myheader;
    myheader = (TMyFileHeader*)((TKey*)ifile_header.GetListOfKeys()->At(0))->ReadObj();
    int tstart = myheader->GetStartTime();
    int tend = myheader->GetEndTime();
    ifile_header.Close();
    
    TFile ifile(string(charging_up_Output+cut_type+"/"+to_string(run)+"/charging_up.root").data(),"READ");
    TGraphErrors *dummygraph;
    TF1 *dummyfunc;
    
    ifile.GetObject("_view0",dummygraph);
    if(!dummygraph){
      cout << "ERROR: no graph named view0 in run " << run << " at cut " << cut_type << endl;
    }
    else if(dummygraph->GetN() > 5){
      dummyfunc = (TF1*)dummygraph->GetListOfFunctions()->At(0);
      slope_0.first.push_back((int)(tend+tstart)/2);
      slope_0.second.push_back(dummyfunc->GetParameter(0));
      intercept_0.first.push_back((int)(tend+tstart)/2);
      intercept_0.second.push_back(dummyfunc->GetParameter(1));
    }
    dummygraph = 0;
    dummyfunc = 0;
    ifile.GetObject("_view1",dummygraph);
    if(!dummygraph){
      cout << "ERROR: no graph named view1 in run " << run << " at cut " << cut_type << endl;
    }
    else if(dummygraph->GetN() > 5){
      dummyfunc = (TF1*)dummygraph->GetListOfFunctions()->At(0);
      slope_1.first.push_back((int)(tend+tstart)/2);
      slope_1.second.push_back(dummyfunc->GetParameter(0));
      intercept_1.first.push_back((int)(tend+tstart)/2);
      intercept_1.second.push_back(dummyfunc->GetParameter(1));
    }
    dummygraph = 0;
    dummyfunc = 0;
    
    for(auto lem : lems){
      ifile.GetObject(string("_LEM_"+to_string(lem)+"_view0").data(),dummygraph);
      if(!dummygraph){
        cout << "ERROR: no graph named " << string("_LEM_"+to_string(lem)+"_view0") << " in run " << run << " at cut " << cut_type << endl;
      }
      else if(dummygraph->GetN() > 5){
        dummyfunc = (TF1*)dummygraph->GetListOfFunctions()->At(0);
        slope_0_ByLEMs[lem].first.push_back((int)(tend+tstart)/2);
        slope_0_ByLEMs[lem].second.push_back(dummyfunc->GetParameter(0));
        intercept_0_ByLEMs[lem].first.push_back((int)(tend+tstart)/2);
        intercept_0_ByLEMs[lem].second.push_back(dummyfunc->GetParameter(1));
      }
      dummygraph = 0;
      dummyfunc = 0;
      ifile.GetObject(string("_LEM_"+to_string(lem)+"_view1").data(),dummygraph);
      if(!dummygraph){
        cout << "ERROR: no graph named " << string("_LEM_"+to_string(lem)+"_view1") << " in run " << run << " at cut " << cut_type << endl;
      }
      else if(dummygraph->GetN() > 5){
        dummyfunc = (TF1*)dummygraph->GetListOfFunctions()->At(0);
        slope_1_ByLEMs[lem].first.push_back((int)(tend+tstart)/2);
        slope_1_ByLEMs[lem].second.push_back(dummyfunc->GetParameter(0));
        intercept_1_ByLEMs[lem].first.push_back((int)(tend+tstart)/2);
        intercept_1_ByLEMs[lem].second.push_back(dummyfunc->GetParameter(1));
      }
      dummygraph = 0;
      dummyfunc = 0;
    }
    ifile.Close();

  }
  
  //slope
  string timeformatshort = "%m/%d";
  string timeformatlong = "#splitline{%m/%d}{%Hh%M}";
  
  TGraph gr_slope_view0(slope_0.first.size(), &slope_0.first[0], &slope_0.second[0]);
  gr_slope_view0.Draw("A*");
  gr_slope_view0.GetXaxis()->SetTimeDisplay(1);
  gr_slope_view0.GetXaxis()->SetNdivisions(-slope_0.first.size()+1);
  if(slope_0.first[1]-slope_0.first[0] > 432000){
    gr_slope_view0.GetXaxis()->SetTimeFormat(timeformatshort.data());
  }
  else{
    gr_slope_view0.GetXaxis()->SetTimeFormat(timeformatlong.data());//"#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00"
  }
  gr_slope_view0.GetXaxis()->SetLabelOffset(.02);
  gr_slope_view0.Draw("A*");
  gPad->Update();
  gPad->Modified();
  gPad->SaveAs(string(charging_up_Output+cut_type+"/slope_view0.png").data());
  
  TGraph gr_slope_view1(slope_1.first.size(), &slope_1.first[0], &slope_1.second[0]);
  gr_slope_view1.Draw("A*");
  gr_slope_view1.GetXaxis()->SetTimeDisplay(1);
  gr_slope_view1.GetXaxis()->SetNdivisions(-slope_1.first.size()+1);
  if(slope_1.first[1]-slope_1.first[0] > 432000){
    gr_slope_view1.GetXaxis()->SetTimeFormat(timeformatshort.data());
  }
  else{
    gr_slope_view1.GetXaxis()->SetTimeFormat(timeformatlong.data());
  }
  gr_slope_view1.GetXaxis()->SetLabelOffset(.02);
  gr_slope_view1.Draw("A*");
  gPad->Update();
  gPad->Modified();
  gPad->SaveAs(string(charging_up_Output+cut_type+"/slope_view1.png").data());
  
  for(auto lem : lems){
    
    TGraph gr_slope_view0_lem(slope_0_ByLEMs[lem].first.size(), &slope_0_ByLEMs[lem].first[0], &slope_0_ByLEMs[lem].second[0]);
    gr_slope_view0_lem.Draw("A*");
    gr_slope_view0_lem.GetXaxis()->SetTimeDisplay(1);
    gr_slope_view0_lem.GetXaxis()->SetNdivisions(-slope_0_ByLEMs[lem].first.size()+1);
    if(slope_0_ByLEMs[lem].first[1]-slope_0_ByLEMs[lem].first[0] > 432000){
      gr_slope_view0_lem.GetXaxis()->SetTimeFormat(timeformatshort.data());
    }
    else{
      gr_slope_view0_lem.GetXaxis()->SetTimeFormat(timeformatlong.data());
    }
    gr_slope_view0_lem.GetXaxis()->SetLabelOffset(.02);
    gr_slope_view0_lem.Draw("A*");
    gPad->Update();
    gPad->Modified();
    gPad->SaveAs(string(charging_up_Output+cut_type+"/slope_view0_LEMS_"+to_string(lem)+".png").data());
    
    TGraph gr_slope_view1_lem(slope_1_ByLEMs[lem].first.size(), &slope_1_ByLEMs[lem].first[0], &slope_1_ByLEMs[lem].second[0]);
    gr_slope_view1_lem.Draw("A*");
    gr_slope_view1_lem.GetXaxis()->SetTimeDisplay(1);
    gr_slope_view1_lem.GetXaxis()->SetNdivisions(-slope_1_ByLEMs[lem].first.size()+1);
    if(slope_1_ByLEMs[lem].first[1]-slope_1_ByLEMs[lem].first[0] > 432000){
      gr_slope_view1_lem.GetXaxis()->SetTimeFormat(timeformatshort.data());
    }
    else{
      gr_slope_view1_lem.GetXaxis()->SetTimeFormat(timeformatlong.data());
    }
    gr_slope_view1_lem.GetXaxis()->SetLabelOffset(.02);
    gr_slope_view1_lem.Draw("A*");
    gPad->Update();
    gPad->Modified();
    gPad->SaveAs(string(charging_up_Output+cut_type+"/slope_view1_LEMS_"+to_string(lem)+".png").data());
  }
  
  
  //intercept
  
  TGraph gr_intercept_view0(intercept_0.first.size(), &intercept_0.first[0], &intercept_0.second[0]);
  gr_intercept_view0.Draw("A*");
  gr_intercept_view0.GetXaxis()->SetTimeDisplay(1);
  gr_intercept_view0.GetXaxis()->SetNdivisions(-intercept_0.first.size()+1);
  if(intercept_0.first[1]-intercept_0.first[0] > 432000){
    gr_intercept_view0.GetXaxis()->SetTimeFormat(timeformatshort.data());
  }
  else{
    gr_intercept_view0.GetXaxis()->SetTimeFormat(timeformatlong.data());
  }
  gr_intercept_view0.GetXaxis()->SetLabelOffset(.02);
  gr_intercept_view0.Draw("A*");
  gPad->Update();
  gPad->Modified();
  gPad->SaveAs(string(charging_up_Output+cut_type+"/intercept_view0.png").data());
  
  TGraph gr_intercept_view1(intercept_1.first.size(), &intercept_1.first[0], &intercept_1.second[0]);
  gr_intercept_view1.Draw("A*");
  gr_intercept_view1.GetXaxis()->SetTimeDisplay(1);
  gr_intercept_view1.GetXaxis()->SetNdivisions(-intercept_1.first.size()+1);
  if(intercept_1.first[1]-intercept_1.first[0] > 432000){
    gr_intercept_view1.GetXaxis()->SetTimeFormat(timeformatshort.data());
  }
  else{
    gr_intercept_view1.GetXaxis()->SetTimeFormat(timeformatlong.data());
  }
  gr_intercept_view1.GetXaxis()->SetLabelOffset(.02);
  gr_intercept_view1.Draw("A*");
  gPad->Update();
  gPad->Modified();
  gPad->SaveAs(string(charging_up_Output+cut_type+"/intercept_view1.png").data());
  
  for(auto lem : lems){
    
    TGraph gr_intercept_view0_lem(intercept_0_ByLEMs[lem].first.size(), &intercept_0_ByLEMs[lem].first[0], &intercept_0_ByLEMs[lem].second[0]);
    gr_intercept_view0_lem.Draw("A*");
    gr_intercept_view0_lem.GetXaxis()->SetTimeDisplay(1);
    gr_intercept_view0_lem.GetXaxis()->SetNdivisions(-intercept_0_ByLEMs[lem].first.size()+1);
    if(intercept_0_ByLEMs[lem].first[1]-intercept_0_ByLEMs[lem].first[0] > 432000){
      gr_intercept_view0_lem.GetXaxis()->SetTimeFormat(timeformatshort.data());
    }
    else{
      gr_intercept_view0_lem.GetXaxis()->SetTimeFormat(timeformatlong.data());
    }
    gr_intercept_view0_lem.GetXaxis()->SetLabelOffset(.02);
    gr_intercept_view0_lem.Draw("A*");
    gPad->Update();
    gPad->Modified();
    gPad->SaveAs(string(charging_up_Output+cut_type+"/intercept_view0_LEMS_"+to_string(lem)+".png").data());
    
    TGraph gr_intercept_view1_lem(intercept_1_ByLEMs[lem].first.size(), &intercept_1_ByLEMs[lem].first[0], &intercept_1_ByLEMs[lem].second[0]);
    gr_intercept_view1_lem.Draw("A*");
    gr_intercept_view1_lem.GetXaxis()->SetTimeDisplay(1);
    gr_intercept_view1_lem.GetXaxis()->SetNdivisions(-intercept_1_ByLEMs[lem].first.size()+1);
    if(intercept_1_ByLEMs[lem].first[1]-intercept_1_ByLEMs[lem].first[0] > 432000){
      gr_intercept_view1_lem.GetXaxis()->SetTimeFormat(timeformatshort.data());
    }
    else{
      gr_intercept_view1_lem.GetXaxis()->SetTimeFormat(timeformatlong.data());
    }
    gr_intercept_view1_lem.GetXaxis()->SetLabelOffset(.02);
    gr_intercept_view1_lem.Draw("A*");
    gPad->Update();
    gPad->Modified();
    gPad->SaveAs(string(charging_up_Output+cut_type+"/intercept_view1_LEMS_"+to_string(lem)+".png").data());
  }
  
  delete gPad;
  
  return;
}
