////////////////////////////////////////////////////////////////////////////////
//
// Produce gain curve using the MPV of the dQdx distribution
// Takes as input the root file produced by utils/track_cuts.cc
// Produce a root file with all the dQdx distributions for every Electric fields
// and a TGraphErrors with the gain curve
//
//TODO: substitute to_string() in histname
//TODO: add scan type selection (if induction, extraction, etc...)
//
// Usage:
// Toggle the scan_type (Extraction, Induction, Amplification) and scan_num (1,2)
// to run the macro and producte the MPV vs field of the scan
//
//
////////////////////////////////////////////////////////////////////////////////

#include <311Lib.h>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <TMyFileHeader.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}



void YZ_gain(int run = 842, string cut_type = "Ds", string v = "July", string m_dQ = "sum", string m_ds = "local"){
  method_ds = m_ds;
  method_dQ = m_dQ;
  string cut_type_and_methods = cut_type + "_" + method_ds + "_" + method_dQ;
  if(!Load_Version(v)){return;}
  gErrorIgnoreLevel = kError;
//  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(0);

//  gStyle->SetPalette(kColorPrintableOnGrey);

  if(!load_cosmics()){return;}
  
  int nbinsY = (int)(2*lem_size/dy);
  int nbinsZ = (int)(6*lem_size/dz);
  vector<int> begins = {-(int)(lem_size/dy), 0};
  vector<int> ends = {(int)(lem_size/dy), (int)(6*lem_size/dz)};
  
  TH2D mpv_crp("mpv","mpv;y;z;gain",nbinsZ,begins[1]*dz,ends[1]*dz,nbinsY,begins[0]*dy,ends[0]*dy);
  TH2D mpv_view0_crp("mpv_view0","mpv_view0;y;z;gain",nbinsZ,begins[1]*dz,ends[1]*dz,nbinsY,begins[0]*dy,ends[0]*dy);
  TH2D mpv_view1_crp("mpv_view1","mpv_view1;y;z;gain",nbinsZ,begins[1]*dz,ends[1]*dz,nbinsY,begins[0]*dy,ends[0]*dy);
  TH1D distri_crp("distri","distri;gain;#",50,0,5);
  TH1D distri_view0_crp("distri_view0","distri_view0;gain;#",50,0,5);
  TH1D distri_view1_crp("distri_view1","distri_view1;gain;#",50,0,5);
  
  map<int,TH2D> mpv; map<int,TH2D> mpv_view0; map<int,TH2D> mpv_view1; map<int,TH1D> distri; map<int,TH1D> distri_view0; map<int,TH1D> distri_view1;
  TH1D stddev("relat_std_dev","relat_std_dev",20,0,15);
  TH1D stddev_view0("relat_std_dev_view0","relat_std_dev_view0",20,0,15);
  TH1D stddev_view1("relat_std_dev_view1","relat_std_dev_view1",20,0,15);
  
  string name, title;
  
  for(auto lem : lems){
    begins = find_YZ(lem);
    begins[0] = begins[0]/dy;
    begins[1] = begins[1]/dz;
    ends[0] = begins[0]+(int)(lem_size/dy);
    ends[1] = begins[1]+(int)(lem_size/dz);
    nbinsY = (int)(lem_size/dy);
    nbinsZ = (int)(lem_size/dz);
    
    name = "mpv_LEM_" + to_string(lem);
    title = "mpv_LEM_" + to_string(lem)+";y;z;gain";
    mpv[lem] = TH2D(name.data(),title.data(),nbinsZ,begins[1]*dz,ends[1]*dz,nbinsY,begins[0]*dy,ends[0]*dy);
    
    name = "distri_LEM_" + to_string(lem);
    title = "distri_LEM_" + to_string(lem)+";gain;#";
    distri[lem] = TH1D(name.data(),title.data(),50,0,5);
    
    name = "mpv_view0_LEM_" + to_string(lem);
    title = "mpv_view0_LEM_" + to_string(lem)+";y;z;gain";
    mpv_view0[lem] = TH2D(name.data(),title.data(),nbinsZ,begins[1]*dz,ends[1]*dz,nbinsY,begins[0]*dy,ends[0]*dy);
    
    name = "distri_view0_LEM_" + to_string(lem);
    title = "distri_view0_LEM_" + to_string(lem)+";gain;#";
    distri_view0[lem] = TH1D(name.data(),title.data(),50,0,5);
    
    name = "mpv_view1_LEM_" + to_string(lem);
    title = "mpv_view1_LEM_" + to_string(lem)+";y;z;gain";
    mpv_view1[lem] = TH2D(name.data(),title.data(),nbinsZ,begins[1]*dz,ends[1]*dz,nbinsY,begins[0]*dy,ends[0]*dy);
    
    name = "distri_view1_LEM_" + to_string(lem);
    title = "distri_view1_LEM_" + to_string(lem)+";gain;#";
    distri_view1[lem] = TH1D(name.data(),title.data(),50,0,5);
  }
  
  nbinsY = (int)(2*lem_size/dy);
  nbinsZ = (int)(6*lem_size/dz);
  begins = {-(int)(lem_size/dy), 0};
  ends = {(int)(lem_size/dy), (int)(6*lem_size/dz)};

  TMyFileHeader header = load_run_header(run);
  if(header.GetRun() == -1){return;}
  if(header.GetYZGainsEff().find(cut_type_and_methods) == header.GetYZGainsEff().end()){
    cout << "ERROR: cut type " << cut_type_and_methods << " not found in run " << run << endl;
    return;
  }
  #if verbose 
  cout << "  Reading run " << run << " header file..." << endl;
  #endif
  auto gains = header.GetYZGainsEff()[cut_type_and_methods];
  auto mpvs = header.GetYZMPVs()[cut_type_and_methods];
  auto processed = header.AreYZGainsProcessed()[cut_type_and_methods];
  
  
  for( int Y = begins[0]; Y < ends[0]; Y++ ){
    double y = (Y+0.5)*dy;
    for( int Z = begins[1]; Z < ends[1]; Z++ ){
      double z = (Z+0.5)*dz;
      pair<int,int> YZ = make_pair(Y,Z);
      if(!processed[YZ]){return;}
      double gain0 = mpvs[YZ].first/mpv_cosmics;
      double gain1 = mpvs[YZ].second/mpv_cosmics;
      double gain_eff = gains[YZ];
      int lem = find_lem(y,z);
      
      mpv_view0_crp.Fill(z, y, gain0);
      distri_view0_crp.Fill(gain0);
      mpv_view0[lem].Fill(z, y, gain0);
      distri_view0[lem].Fill(gain0);
    
      mpv_view1_crp.Fill(z, y, gain1);
      distri_view1_crp.Fill(gain1);
      mpv_view1[lem].Fill(z, y, gain1);
      distri_view1[lem].Fill(gain1);
    
      mpv_crp.Fill(z, y, gain_eff);
      distri_crp.Fill(gain_eff);
      mpv[lem].Fill(z, y, gain_eff);
      distri[lem].Fill(gain_eff);
      
    }//for z
  }// for y
  #if verbose
  cout << "    done reading histograms and filling graphs" << endl;
  #endif
  
  string outpath, outfile;
  int cany = 850;
  int canz = 1000;
  
  for(auto lem : lems){
    nbinsY = (int)(lem_size/dy);
    nbinsZ = (int)(lem_size/dz);
    outpath = YZ_Output;
    check_and_mkdir(outpath);
    outpath = outpath + cut_type;
    check_and_mkdir(outpath);
    outpath = outpath + "/" + to_string(run);
    check_and_mkdir(outpath);
    outpath = outpath + "/" + to_string(lem) + "/";
    check_and_mkdir(outpath);
    outfile = outpath + "gain.root";
    TFile ofilelem(outfile.data(), "RECREATE");
    #if verbose
    cout << "    Writing file: " << outfile << endl;
    #endif
    
    
    TCanvas *can1 = new TCanvas(string(string(mpv[lem].GetName())+"_pad").data(),string(string(mpv[lem].GetName())+"_pad").data(),canz,cany);
    gStyle->SetOptStat(0);
    ofilelem.cd();
    mpv[lem].SetMinimum(0);
    mpv[lem].SetMaximum(5);
    mpv[lem].GetXaxis()->SetNdivisions(-nbinsZ);
    mpv[lem].GetXaxis()->SetLabelSize(0.02);
    mpv[lem].GetYaxis()->SetNdivisions(-nbinsY);
    mpv[lem].Draw("COLZ");
    mpv[lem].Write();
    can1->Write();
    can1->SaveAs(string(outpath + string(mpv[lem].GetName()) + ".png").data());
    delete can1;
    
    TCanvas *can2 = new TCanvas(string(string(mpv_view0[lem].GetName())+"_pad").data(),string(string(mpv_view0[lem].GetName())+"_pad").data(),canz,cany);
    ofilelem.cd();
    mpv_view0[lem].SetMinimum(0);
    mpv_view0[lem].SetMaximum(5);
    mpv_view0[lem].GetXaxis()->SetNdivisions(-nbinsZ);
    mpv_view0[lem].GetXaxis()->SetLabelSize(0.02);
    mpv_view0[lem].GetYaxis()->SetNdivisions(-nbinsY);
    mpv_view0[lem].Draw("COLZ");
    mpv_view0[lem].Write();
    can2->Write();
    can2->SaveAs(string(outpath + string(mpv_view0[lem].GetName()) + ".png").data());
    delete can2;
    
    TCanvas *can3 = new TCanvas(string(string(mpv_view1[lem].GetName())+"_pad").data(),string(string(mpv_view1[lem].GetName())+"_pad").data(),canz,cany);
    ofilelem.cd();
    mpv_view1[lem].SetMinimum(0);
    mpv_view1[lem].SetMaximum(5);
    mpv_view1[lem].GetXaxis()->SetNdivisions(-nbinsZ);
    mpv_view1[lem].GetXaxis()->SetLabelSize(0.02);
    mpv_view1[lem].GetYaxis()->SetNdivisions(-nbinsY);
    mpv_view1[lem].Draw("COLZ");
    mpv_view1[lem].Write();
    can3->Write();
    can3->SaveAs(string(outpath + string(mpv_view1[lem].GetName()) + ".png").data());
    delete can3;
    
    gStyle->SetOptStat(1111);
    ofilelem.cd();
    distri[lem].Draw();
    distri[lem].Write();
    gPad->SetName(string(string(distri[lem].GetName())+"_pad").data());
    gPad->Write();
    gPad->SaveAs(string(outpath + string(distri[lem].GetName()) + ".png").data());
    stddev.Fill(100*distri[lem].GetStdDev()/distri[lem].GetMean());
    
    ofilelem.cd();
    distri_view0[lem].Draw();
    distri_view0[lem].Write();
    gPad->SetName(string(string(distri_view0[lem].GetName())+"_pad").data());
    gPad->Write();
    gPad->SaveAs(string(outpath + string(distri_view0[lem].GetName()) + ".png").data());
    stddev_view0.Fill(100*distri_view0[lem].GetStdDev()/distri[lem].GetMean());
    
    ofilelem.cd();
    distri_view1[lem].Draw();
    distri_view1[lem].Write();
    gPad->SetName(string(string(distri_view1[lem].GetName())+"_pad").data());
    gPad->Write();
    gPad->SaveAs(string(outpath + string(distri_view1[lem].GetName()) + ".png").data());
    stddev_view1.Fill(100*distri_view1[lem].GetStdDev()/distri[lem].GetMean());
    delete gPad;
    
    ofilelem.Close();
  }
  
  nbinsY = (int)(2*lem_size/dy);
  nbinsZ = (int)(6*lem_size/dz);
  outpath = YZ_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type;
  check_and_mkdir(outpath);
  outpath = outpath + "/" + to_string(run);
  check_and_mkdir(outpath);
  outpath = outpath + "/crp/";
  check_and_mkdir(outpath);
  outfile = outpath + "gain.root";
  TFile ofile(outfile.data(), "RECREATE");
  #if verbose
  cout << "    Writing file: " << outfile << endl;
  #endif
  
  canz = 3000;
  vector<TLine*> lines = {};
  lines.push_back(new TLine(0,0,6*48,0));
  lines.push_back(new TLine(48,-48,48,48));
  lines.push_back(new TLine(2*48,-48,2*48,48));
  lines.push_back(new TLine(3*48,-48,3*48,48));
  lines.push_back(new TLine(4*48,-48,4*48,48));
  lines.push_back(new TLine(5*48,-48,5*48,48));
  
  TCanvas *can1 = new TCanvas(string(string(mpv_crp.GetName())+"_pad").data(),string(string(mpv_crp.GetName())+"_pad").data(),canz,cany);
  gStyle->SetOptStat(0);
  ofile.cd();
  mpv_crp.SetMinimum(0);
  mpv_crp.SetMaximum(5);
  mpv_crp.GetXaxis()->SetNdivisions(-nbinsZ);
  mpv_crp.GetXaxis()->SetLabelSize(0.02);
  mpv_crp.GetYaxis()->SetNdivisions(-nbinsY);
  mpv_crp.Draw("COLZ");
  for(auto myline : lines){
    myline->SetLineWidth(2);
    myline->Draw();
  }
  mpv_crp.Write();
  can1->Write();
  can1->SaveAs(string(outpath + string(mpv_crp.GetName()) + ".png").data());
  delete can1;
  
  TCanvas *can2 = new TCanvas(string(string(mpv_view0_crp.GetName())+"_pad").data(),string(string(mpv_view0_crp.GetName())+"_pad").data(),canz,cany);
  ofile.cd();
  mpv_view0_crp.SetMinimum(0);
  mpv_view0_crp.SetMaximum(5);
  mpv_view0_crp.GetXaxis()->SetNdivisions(-nbinsZ);
  mpv_view0_crp.GetXaxis()->SetLabelSize(0.02);
  mpv_view0_crp.GetYaxis()->SetNdivisions(-nbinsY);
  mpv_view0_crp.Draw("COLZ");
  for(auto myline : lines){
    myline->Draw();
  }
  mpv_view0_crp.Write();
  can2->Write();
  can2->SaveAs(string(outpath + string(mpv_view0_crp.GetName()) + ".png").data());
  delete can2;
  
  TCanvas *can3 = new TCanvas(string(string(mpv_view1_crp.GetName())+"_pad").data(),string(string(mpv_view1_crp.GetName())+"_pad").data(),canz,cany);
  ofile.cd();
  mpv_view1_crp.SetMinimum(0);
  mpv_view1_crp.SetMaximum(5);
  mpv_view1_crp.GetXaxis()->SetNdivisions(-nbinsZ);
  mpv_view1_crp.GetXaxis()->SetLabelSize(0.02);
  mpv_view1_crp.GetYaxis()->SetNdivisions(-nbinsY);
  mpv_view1_crp.Draw("COLZ");
  for(auto myline : lines){
    myline->Draw();
  }
  mpv_view1_crp.Write();
  can3->Write();
  can3->SaveAs(string(outpath + string(mpv_view1_crp.GetName()) + ".png").data());
  delete can3;
  
  gStyle->SetOptStat(1111);
  ofile.cd();
  distri_crp.Draw();
  distri_crp.Write();
  gPad->SetName(string(string(distri_crp.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(distri_crp.GetName()) + ".png").data());
  
  ofile.cd();
  distri_view0_crp.Draw();
  distri_view0_crp.Write();
  gPad->SetName(string(string(distri_view0_crp.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(distri_view0_crp.GetName()) + ".png").data());
  
  ofile.cd();
  distri_view1_crp.Draw();
  distri_view1_crp.Write();
  gPad->SetName(string(string(distri_view1_crp.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(distri_view1_crp.GetName()) + ".png").data());

  ofile.cd();
  stddev.Draw();
  stddev.Write();
  gPad->SetName(string(string(stddev.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(stddev.GetName()) + ".png").data());
  
  ofile.cd();
  stddev_view0.Draw();
  stddev_view0.Write();
  gPad->SetName(string(string(stddev_view0.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(stddev_view0.GetName()) + ".png").data());
  
  ofile.cd();
  stddev_view1.Draw();
  stddev_view1.Write();
  gPad->SetName(string(string(stddev_view0.GetName())+"_pad").data());
  gPad->Write();
  gPad->SaveAs(string(outpath + string(stddev_view1.GetName()) + ".png").data());
  
  delete gPad;
  for(auto line : lines){delete line;}
  
  ofile.Close();
  
  #if verbose
  cout << "all done!" << endl;
  #endif
  return;
}//end macro
