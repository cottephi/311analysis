#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <TMyFileHeader.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}

vector<string> params = {"date","TE0037","TE0038","TE0039","TE0040","TE0041","TE0042","TE0043","TE0044","TE0045","TE0046","TE0047","TE0048","TE1001","TE1002","TE1003","TE1004","TE0049","TE0050","TE0051","TE0052","TE0053","TE0054","TE0055","TE0056","TE0057","TE0058","TE0059","TE0060","PE0006"};
//vector<string> params = {"date","TE0037","TE0041","TE0045","TE1001","PE0006"};

void pressure(string srun = "801"){
  
  gErrorIgnoreLevel = kError;

  vector <string> temp_crp = {"TE0037","TE0038","TE0039","TE0040","TE0041","TE0042","TE0043","TE0044","TE0045","TE0046","TE0047","TE0048","TE1001","TE1002","TE1003","TE1004"};
  map<string,double> temp_cryostat;
  temp_cryostat["TE0049"] = 115;
  temp_cryostat["TE0050"] = 119;
  temp_cryostat["TE0051"] = 123;
  temp_cryostat["TE0052"] = 127;
  temp_cryostat["TE0053"] = 131;
  temp_cryostat["TE0054"] = 135;
  temp_cryostat["TE0055"] = 139;
  temp_cryostat["TE0056"] = 143;
  temp_cryostat["TE0057"] = 147;
  temp_cryostat["TE0058"] = 151;
  temp_cryostat["TE0059"] = 155;
  temp_cryostat["TE0060"] = 159;
  string directory = slow_control+srun+"/";
  string path = directory;
  if(srun.find("all") == string::npos){path = path+srun+".txt";}
  else{path = path+srun+".dat";}
  if(!ExistTest(path)){cout << "ERROR: file " << path << " not found" << endl; return;}
  if(!ExistTest(string(directory+"pressure_and_temp/").data())){
    #if verbose
    cout << "Creating directory " << directory << "pressure_and_temp..." << endl;
    #endif
    mkdir(string(directory+"pressure_and_temp").data(),0777);
  }
  directory = directory + "pressure_and_temp/";
  string command = "sed 's/[[:blank:]]/;/g' " + path + " > tmp.txt";
  system(command.data());
  command = "perl -p -e 's/;\\r\\n/\\n/' tmp.txt > tmp2.txt";
  system(command.data());

  vector<string> params_in_file = {};
  fstream ifile;
  ifile.open("tmp2.txt", fstream::in);
  if(!ifile.is_open()){
    cout << "ERROR : can not open file tmp2.txt" << endl;
    return;
  }
  int i = 0;
  vector<int> columns_to_rec = {};
  map<string,vector<double> > map_par_values;
  vector<double> mean_cryo_temps = {};
  vector<double> cryo_pos_for_plot = {};
  double rho_var = 0;
  double p_mean = 0;
  
  while(!ifile.eof()){
    char buffer[1000];
    ifile.getline(buffer, 1000);
    char* tokBuffer;
    if( string(buffer).find(";") == string::npos ){break;}
    tokBuffer = strtok(buffer, ";");
    int col = 0;
    if(i == 0){
      if(std::find(params.begin(), params.end(), string(tokBuffer)) != params.end()){
        #if verbose
        cout << "reading parameter " << tokBuffer << endl;
        #endif
        params_in_file.push_back(string(tokBuffer));
        columns_to_rec.push_back(col);
        map_par_values[string(tokBuffer)] = {};
      }
      tokBuffer = strtok(NULL, ";");
      while(tokBuffer){
        col++;
        if(std::find(params.begin(), params.end(), string(tokBuffer)) != params.end()){
          #if verbose
          cout << "reading parameter " << tokBuffer << endl;
          #endif
          params_in_file.push_back(string(tokBuffer));
          columns_to_rec.push_back(col);
          map_par_values[string(tokBuffer)] = {};
        }
        tokBuffer = strtok(NULL, ";");
      }
      i++;
      continue;
    }
    if(std::find(columns_to_rec.begin(), columns_to_rec.end(), col) != columns_to_rec.end()){
      map_par_values[params_in_file[std::find(columns_to_rec.begin(), columns_to_rec.end(), col)-columns_to_rec.begin()]].push_back(atof(tokBuffer));
    }
    tokBuffer = strtok(NULL, ";");
    while(tokBuffer){
      col++;
      if(std::find(columns_to_rec.begin(), columns_to_rec.end(), col) != columns_to_rec.end()){
        map_par_values[params_in_file[std::find(columns_to_rec.begin(), columns_to_rec.end(), col)-columns_to_rec.begin()]].push_back(atof(tokBuffer));
      }
      tokBuffer = strtok(NULL, ";");
    }
  }
  ifile.close();
  
  if(map_par_values["data"].size() == 0){
    cout << "Error in pressure.cc: empty run" << endl;
    map_par_values.clear();
    return;
  }
  
  cout << "...done" << endl;
  TFile ofile_pressure(string(directory + "pressure.root").data(),"RECREATE");
  TFile ofile_temperature_CRP(string(directory + "temperature_CRP.root").data(),"RECREATE");
  TFile ofile_temperature_cryostat(string(directory + "temperature_cryostat.root").data(),"RECREATE");
  
  
  for(auto par : map_par_values){
    int i=0;
////////////////////////////////////////////////
    //param vs time
////////////////////////////////////////////////
    cout << "Plotting " << par.first << "..." << endl;

    if(par.first == "date"){continue;}
    i++;
    cout << i << endl;
    TGraph gr(par.second.size(), &map_par_values["date"][0], &par.second[0]);
    i++;
    cout << i << endl;
    gr.SetTitle(string(par.first).data());
    i++;
    cout << i << endl;
    gr.SetName(string("graph_"+par.first).data());
    i++;
    cout << i << endl;
    gr.Draw("AC");
    i++;
    cout << i << endl;
    gr.GetXaxis()->SetTimeDisplay(1);
    i++;
    cout << i << endl;
    gr.GetXaxis()->LabelsOption("v");
    i++;
    cout << i << endl;
    if(map_par_values["date"].back()-map_par_values["date"][0] > 24*3600){
    i++;
    cout << " " << i << endl;
      gr.GetXaxis()->SetTimeFormat("%d/%m %Hh%M");
      gr.GetXaxis()->SetNdivisions(5);
    }
    else{
    i++;
    cout << "  " << i << endl;
      gr.GetXaxis()->SetTimeFormat("%Hh%M");
      gr.GetXaxis()->SetNdivisions(8);
    }
    i++;
    cout << i << endl;
    gr.Draw("AC");
    i++;
    cout << i << endl;
    gPad->Update();
    i++;
    cout << i << endl;
    gPad->Modified();
    i++;
    cout << i << endl;
    if(par.first == "PE0006"){ofile_pressure.cd();}
    else if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
    else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
    i++;
    cout << i << endl;
    gr.Write();
    i++;
    cout << i << endl;
    gPad->SetName(string(string(gr.GetName()) + "_pad").data());
    i++;
    cout << i << endl;
//    gPad->SaveAs(string(directory + "graph_"+par.first + ".png").data());
    if(par.first == "PE0006"){ofile_pressure.cd();}
    else if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
    else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
    gPad->Write();
    delete gPad;
    
////////////////////////////////////////////////
    //param distribution
////////////////////////////////////////////////

    TH1D h(string("histo_"+par.first).data(),string(par.first+";"+par.first+";#").data(),100,0.99*par.second[min_element(par.second.begin(),par.second.end()) - par.second.begin()],par.second[max_element(par.second.begin(),par.second.end()) - par.second.begin()]*1.01);
    for(auto value : par.second){
      h.Fill(value);
    }
    if(temp_cryostat.find(par.first) != temp_cryostat.end()){cryo_pos_for_plot.push_back(temp_cryostat[par.first]); mean_cryo_temps.push_back(h.GetMean());}
    h.Draw();
    if(par.first == "PE0006"){ofile_pressure.cd();}
    else if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
    else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
    h.Write();
    gPad->SetName(string(string(h.GetName()) + "_pad").data());
//    gPad->SaveAs(string(directory + "histo_"+par.first + ".png").data());
    if(par.first == "PE0006"){ofile_pressure.cd();p_mean = h.GetMean();}
    else if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
    else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
    gPad->Write();
    delete gPad;
    
    if(par.first != "PE0006"){
    
////////////////////////////////////////////////
    //density vs time
////////////////////////////////////////////////

      vector<double> rho = {};
      for(int i = 0; i < par.second.size(); i++){
        rho.push_back(map_par_values["PE0006"][i]/par.second[i]);
      }
      TGraph gr_rho(rho.size(), &map_par_values["date"][0], &rho[0]);
      gr_rho.SetTitle(string("PE0006/"+par.first).data());
      gr_rho.SetName(string("graph_PE0006Over"+par.first).data());
      gr_rho.Draw("AC");
      if(par.first == "TE0056"){rho_var = *max_element(begin(rho),end(rho))-*min_element(begin(rho),end(rho));}
      gr_rho.GetXaxis()->SetTimeDisplay(1);
      gr_rho.GetXaxis()->LabelsOption("v");
      if(map_par_values["date"].back()-map_par_values["date"][0] > 24*3600){
        gr_rho.GetXaxis()->SetTimeFormat("%d/%m %Hh%M");
        gr_rho.GetXaxis()->SetNdivisions(5);
      }
      else{
        gr_rho.GetXaxis()->SetTimeFormat("%Hh%M");
        gr_rho.GetXaxis()->SetNdivisions(8);
      }
      gr_rho.Draw("AC");
      gPad->Update();
      gPad->Modified();
      if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
      else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
      gr_rho.Write();
      gPad->SetName(string(string(gr_rho.GetName()) + "_pad").data());
//      gPad->SaveAs(string(directory + "graph_PE0006Over"+par.first + ".png").data());
      if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
      else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
      gPad->Write();
      delete gPad;
      
////////////////////////////////////////////////
    //density distribution
////////////////////////////////////////////////
      
      TH1D h_rho(string("histo_PE0006Over"+par.first).data(),string("PE0006Over"+par.first+";PE0006Over"+par.first+";#").data(),100,0.99*rho[min_element(rho.begin(),rho.end()) - rho.begin()],rho[max_element(rho.begin(),rho.end()) - rho.begin()]*1.01);
      for(auto value : rho){
        h_rho.Fill(value);
      }
      h_rho.Draw();
      if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
      else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
      h_rho.Write();
      gPad->SetName(string(string(h_rho.GetName()) + "_pad").data());
//      gPad->SaveAs(string(directory + "histo_PE0006Over"+par.first + ".png").data());
      gPad->Write();
      delete gPad;
    }
    
    cout << "chien" << endl;
  }
  ofile_pressure.Close();
  ofile_temperature_CRP.Close();
      
////////////////////////////////////////////////
    //Computing LEM temperature
////////////////////////////////////////////////
  
  TGraph temp_vs_height(cryo_pos_for_plot.size(), &cryo_pos_for_plot[0], &mean_cryo_temps[0]);
  temp_vs_height.SetTitle("temperature vs height");
  temp_vs_height.SetName("temperature_vs_height");
  temp_vs_height.Draw("A*");
  
  vector<TLatex*> annotations(temp_vs_height.GetN());
  for(int i = 0; i < temp_vs_height.GetN(); i++){
    double x,y;
    string myparam;
    temp_vs_height.GetPoint(i,x,y);
    for(auto item : temp_cryostat){if(item.second == x){myparam=item.first; break;}}
    annotations[i] = new TLatex(temp_vs_height.GetX()[i], temp_vs_height.GetY()[i],myparam.data());
    annotations[i]->SetTextSize(0.02); annotations[i]->SetTextColor(kRed);
    annotations[i]->Draw();
  }
  TF1 *temp_const = new TF1("temp_const","pol0",110,143);
  TF1 *temp_gradient = new TF1("temp_gradient","pol1",143,160);
  temp_vs_height.Fit(temp_gradient,"R+");
  temp_vs_height.Fit(temp_const,"R+");
  double CRP_temp = temp_const->GetParameter(0)+0.5*temp_gradient->GetParameter(1);
  cout << p_mean/CRP_temp << " " << rho_var << endl;
  ofile_temperature_cryostat.cd();
  temp_vs_height.Write();
  gPad->SetName(string(string(temp_vs_height.GetName()) + "_pad").data());
  gPad->SaveAs(string(directory + "graph_"+string(temp_vs_height.GetName()).data() + ".png").data());
  ofile_temperature_cryostat.cd();
  gPad->Write();
  delete gPad;
  
////////////////////////////////////////////////
    //gain variations for this run
////////////////////////////////////////////////

  vector<double> fields_to_plot = {29,30,31,32,33,34,35};
  vector<double> gains_middle;
  vector<double> gains_up;
  vector<double> gains_down;
  for(auto fiel : fields_to_plot){
    gains_middle.push_back(theoretical_gain(1,p_mean/CRP_temp,fiel));
    gains_up.push_back(theoretical_gain(1,p_mean/CRP_temp + rho_var/2.,fiel));
    gains_down.push_back(theoretical_gain(1,p_mean/CRP_temp - rho_var/2.,fiel));
  }
  TGraph gain_middle(gains_middle.size(),&fields_to_plot[0],&gains_middle[0]);
  gain_middle.SetMarkerColor(kGreen);
  gain_middle.SetLineColor(kGreen);
  TGraph gain_up(gains_up.size(),&fields_to_plot[0],&gains_up[0]);
  gain_up.SetMarkerColor(kRed);
  gain_up.SetLineColor(kRed);
  TGraph gain_down(gains_down.size(),&fields_to_plot[0],&gains_down[0]);
  gain_down.SetMarkerColor(kRed);
  gain_down.SetLineColor(kRed);
  gain_up.SetTitle(string("gain with density variating from " + to_string((int)p_mean/CRP_temp - rho_var/2.) + " to " + to_string((int)p_mean/CRP_temp + rho_var/2.)+";LEM Field (kV/cm);Gain in LEM").data());
  gain_up.Draw("AC");
  gain_middle.Draw("SAME C");
  gain_down.Draw("SAME C");
  gPad->SetName("gain_variations");
  gPad->SaveAs(string(directory + "gain_variations.png").data());
  ofile_temperature_cryostat.cd();
  gPad->Write();
  delete gPad;
      
  ofile_temperature_cryostat.Close();
////////////////////////////////////////////////
    //Updating header file
////////////////////////////////////////////////
  
  system("rm ./tmp*");
  if(srun.find("all") != string::npos){
    map_par_values;
    return;
  }
  
  TFile headerfile(string(runs_headers+srun+".root").data(),"UPDATE");
  TMyFileHeader *header = (TMyFileHeader*)((TKey*)headerfile.GetListOfKeys()->At(0))->ReadObj();
  header->SetRho(p_mean/CRP_temp);
  header->SetRhoVar(rho_var);
  header->Write(string("header_"+srun).data(),TObject::kOverwrite);
  headerfile.Close();
  
  return;
  map_par_values;
}
