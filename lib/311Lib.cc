////////////////////////////////////////////////////////////////////////////////
//
// 311 common analysis function and variables
//
////////////////////////////////////////////////////////////////////////////////

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <glob.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TString.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveStats.h>
#include <THStack.h>
#include <TList.h>
#include <TKey.h>
#include <TMyFileHeader.h>
#include <TLatex.h>
//project libraries
#include "/eos/user/p/pcotte/311analysis/lib/311Lib.h"

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}

//Utilities ********************************************************************

bool Load_Version(string v){
  #if verbose
  cout << "Loading version " << v << endl;
  #endif
  version = v;
  if(version == "Feb"){
    path_311data = "/eos/user/p/pcotte/311data/2018_Feb_05/";
    path_wa105_311data = "/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_Feb_05/ROOT/recofast/";
    recofile_suffix = "-Parser.root";
    tpc_boundaries = {-50, 50, -50, 50, 0, 300};
    pitch = 0.3125;
    dx=10;
    dy=5;
    dz=5;
    lem_size = 50;
  }
  else if(version == "June"){
    path_311data = "/eos/user/p/pcotte/311data/2018_June_24/";
    path_wa105_311data = "/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofast/";
    recofile_suffix = "-RecoFast-Parser.root";
    tpc_boundaries = {-50, 50, -48, 48, 0, 288};
    pitch = 0.3;
    dx=10;
    dy=4.8;
    dz=4.8;
    lem_size = 48;
  }
  else if(version == "Junec"){
    path_311data = "/eos/user/p/pcotte/311data/2018_June_24c/";
    path_wa105_311data = "/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24c/ROOT/recofast/";
    recofile_suffix = "-RecoFast-Parser.root";
    tpc_boundaries = {-50, 50, -48, 48, 0, 288};
    pitch = 0.3;
    dx=10;
    dy=4.8;
    dz=4.8;
    lem_size = 48;
  }
  else if(version == "July"){
    path_311data = "/eos/user/p/pcotte/311data/2018_July_30/";
    path_wa105_311data = "/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_July_30/ROOT/recofast/";
    recofile_suffix = "-RecoFast-Parser.root";
    tpc_boundaries = {-50, 50, -50, 50, 0, 300};
    pitch = 0.3125;
    dx=10;
    dy=5;
    dz=5;
    lem_size = 50;
  }
  else if(version == "LArSoft_old"){
    path_311data = "/eos/user/p/pcotte/311data/LArSoft/old/";
    path_wa105_311data = "/eos/user/p/pcotte/311analysis/LArSoft/old/";
    recofile_suffix = "-RecoFast-Parser.root";
    tpc_boundaries = {-50, 50, -50, 50, 0, 300};
    pitch = 0.3125;
    dx=10;
    dy=5;
    dz=5;
    lem_size = 50;
  }
  else if(version == "LArSoft_new"){
    path_311data = "/eos/user/p/pcotte/311data/LArSoft/new/";
    path_wa105_311data = "/eos/user/p/pcotte/311analysis/LArSoft/new/";
    recofile_suffix = "-RecoFast-Parser.root";
    tpc_boundaries = {-50, 50, -50, 50, 0, 300};
    pitch = 0.3125;
    dx=10;
    dy=5;
    dz=5;
    lem_size = 50;
  }
  else{
    cout << "Load version ERROR: Unknown reco version " << version << endl;
    return false;
  }
  MinX = tpc_boundaries[0] + vol_cut[0];
  MaxX = tpc_boundaries[1] - vol_cut[1];
  MinY = tpc_boundaries[2] + vol_cut[2];
  MaxY = tpc_boundaries[3] - vol_cut[3];
  MinZ = tpc_boundaries[4] + vol_cut[4];
  MaxZ = tpc_boundaries[5] - vol_cut[5];
  check_and_mkdir(path_311data);
  SelectTrack_Output = path_311data + "selected_tracks_" + method_ds + "_" + method_dQ + "/";
  SelectTrack_MC_Output = path_MC_output + "selected_tracks_" + method_ds + "_" + method_dQ + "/";
  YZ_SelectTrack_Output = path_311data + "YZ_selected_tracks_" + method_ds + "_" + method_dQ + "/";
  GainAna_Output = path_311data + "gain_ana_" + method_ds + "_" + method_dQ + "/";
  GainAna_dQds_Output = path_311data + "gain_ana_dQds_" + method_ds + "_" + method_dQ + "/";
  Purity_Output = path_311data + "purity_" + method_ds + "_" + method_dQ + "/";
  cuts_analysis_Output = path_311data + "cuts_analysis_" + method_ds + "_" + method_dQ + "/";
  MPV_vs_stuff_Output = path_311data + "MPV_vs_stuff_" + method_ds + "_" + method_dQ + "/";
  MPV_vs_stuff_tracks_Output = path_311data + "MPV_vs_stuff_tracks_" + method_ds + "_" + method_dQ + "/";
  gain_by_lem_Output = path_311data + "gain_by_lem_" + method_ds + "_" + method_dQ + "/";
  YZ_Output = path_311data + "YZ_" + method_ds + "_" + method_dQ + "/";
  dQds_MC_Output = path_MC_output + "dQds/";
  charging_up_Output = path_311data + "charging_up_" + method_ds + "_" + method_dQ + "/";
  gain_stability_Output = path_311data + "gain_stability_" + method_ds + "_" + method_dQ + "/";
  dQds_Output = path_311data + "dQds_" + method_ds + "_" + method_dQ + "/";
  dQds_charging_up_Output = path_311data + "dQds_charging_up_" + method_ds + "_" + method_dQ + "/";
  dQds_YZ_Output = path_311data + "dQds_YZ_" + method_ds + "_" + method_dQ + "/"; 
  
  #if verbose
  cout << "...done." << endl;
  #endif
  return true;
}

///////////////////////////////////////////////////////////////////////////////

bool load_cosmics(){
  TH1D *h_cosmics = new TH1D("dqds_cosmics","dqds_cosmics",100,0,50);
  if( !ExistTest(dqds_cosmics) ){
    cout<< "ERROR: cosmic file not found" << endl;
    return false;
  }
  fstream cosmics_file;
  cosmics_file.open(dqds_cosmics.data(), fstream::in);
  if(!cosmics_file.is_open()){
    cout << "ERROR: could not open cosmics simulation file" << endl;
    return false;
  }
  while(!cosmics_file.eof()){
    float dqdx = 0;
    cosmics_file >> dqdx;
    h_cosmics->Fill(dqdx);
  }
  cosmics_file.close();
  mpv_cosmics = h_cosmics->GetBinCenter(h_cosmics->GetMaximumBin());
  TFile *cosmics_TFile = TFile::Open("/eos/user/p/pcotte/311analysis/cosmics.root", "RECREATE");
  h_cosmics->Write();
  cosmics_TFile->Close();
  delete h_cosmics;
  #if verbose
  cout << "Expected MPV of cosmics muons before drift: " << mpv_cosmics << endl;
  #endif
  return true;
}

///////////////////////////////////////////////////////////////////////////////

bool load_runs_2D(string scan_type, int scan_num, vector<int> &run_list, vector<pair<float,float> > &field, vector<string> &const_fields_names, vector<float> &const_fields_values){

  //Amplification_Extraction
  if( scan_type == "Amplification_Extraction" and scan_num == 1 ){
    // Drift: 330  //induction: 1
    run_list = {741, 743, 744, 745, 747, 766, 767, 768, 769, 770, 774, 776, 778, 779, 780, 781, 782, 1177};
  }
  else if( scan_type == "Amplification_Extraction" and scan_num == 2 ){
    // Drift: 500  //induction: 1
    run_list = {783, 784, 785, 786, 787, 788, 981, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 993, 994, 995, 996, 997, 998, 999, 1000, 1002, 1003, 1004, 1005, 1006, 1199};
  }
  else if( scan_type == "Amplification_Extraction" and scan_num == 3 ){
    // Drift: 500  //induction: 1.5
    run_list = {1198, 1194, 833, 834, 835, 836, 837, 838, 840, 841, 842, 843};
  }
  else if( scan_type == "Amplification_Extraction" and scan_num == 4 ){
    // Drift: 500  //induction: 1.25
    run_list = {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1016, 1035, 1036, 1037, 1038, 1039, 1040};
  }
  else if( scan_type == "Amplification_Extraction" and scan_num == 5 ){
    // Drift: 300  //induction: 2.5
    run_list = {1165, /*1166,*/ 1178, 1180, 1181, 1182, 1183, 1187, 1188};
  }
  else if( scan_type == "Amplification_Extraction" and scan_num == 6 ){
    // Drift: 500  //induction: 2.5
    run_list = {1189, 1190, 1191, 1192, 1195};
  }
  
  //Amplification_Induction
  else if( scan_type == "Amplification_Induction" and scan_num == 1 ){
    // Drift: 500  //Extraction: 2.2
    run_list = {792, 803, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 993, 994, 995, 996, 997, 998, 999, 1000};
  }
  else if( scan_type == "Amplification_Induction" and scan_num == 2 ){
    // Drift: 500  //Extraction: 2.0
    run_list = {790, 801, 1012, 1013, 1014, 1016, 1035, 1036, 1037, 1038, 837, 838};
  }
  else if( scan_type == "Amplification_Induction" and scan_num == 3 ){
    // Drift: 500  //Extraction: 1.9
    run_list = {789, 800, 836, 840, 841, 842, 843};
  }
  else if( scan_type == "Amplification_Induction" and scan_num == 4 ){
    // Drift: 330  //Extraction: 1.2
    run_list = {744, 745, 774, 776, 778, 1178};
  }
  
  //Induction_Extraction
  else if( scan_type == "Induction_Extraction" and scan_num == 1 ){
    // Drift: 500  //Amplification: 28
    run_list = {783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 833, 834, 835, 836, 837, 838, 840, 841, 842, 843, 1039, 1040, 1195};
  }
  else if( scan_type == "Induction_Extraction" and scan_num == 2 ){
    // Drift: 300  //Amplification: 28
    run_list = {770, 774, 776, 778, 779, 780, 781, 782, 765, 769, 747};
  }
  else if( scan_type == "Induction_Extraction" and scan_num == 3 ){
    // Drift: 500  //Amplification: 27
    run_list = {1000, 1002, 1003, 1004, 1007, 1008, 1009, 1010, 1011, 1191};
  }
  else if( scan_type == "Induction_Extraction" and scan_num == 4 ){
    // Drift: 500  //Amplification: 27.5
    run_list = {1005, 1006, 1012, 1013, 1014, 1016, 1035, 1036, 1037, 1038};
  }
  else if( scan_type == "Induction_Extraction" and scan_num == 5 ){
    // Drift: 300  //Amplification: 23
    run_list = {1172, 1173, 1174, 1175, 1176, 1177, 1178, 1180, 1181};
  }
  else if( scan_type == "Induction_Extraction" and scan_num == 6 ){
    // Drift: 500  //Amplification: 25
    run_list = {981, 982, 983, 984, 985, 986, 987, 988, 989, 1190};
  }

  else{
    #if verbose
    cout<< "Error: unknown combination of scan type and scan number (" << scan_type << ", " << scan_num << ")." << endl;
    #endif
    return false;
  }
  string field1 = scan_type.substr(0,scan_type.find_first_of("_"));
  string field2 = scan_type.substr(scan_type.find_first_of("_")+1);
  for(auto r : run_list){
    field.push_back(make_pair(runs_and_fields[r][field1], runs_and_fields[r][field2]));
  }
  for(auto f : runs_and_fields[run_list[0]]){
    if(f.first == field1 or f.first == field2){continue;}
    const_fields_names.push_back(f.first);
    const_fields_values.push_back(f.second);
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool load_runs(string scan_type, int scan_num, vector<int> &run_list, vector<float> &field, vector<string> &const_fields_names, vector<float> &const_fields_values){

  string field_to_read = scan_type;
  if(scan_type == "All"){
    field_to_read = "Amplification";
  }
  else if(scan_type == "Identical"){
    field_to_read = "None";
  }
  
  //Generated with python macro (all runs, 0% drift tolerence)
  // Drift :
  if( scan_type == "Drift" and scan_num == 1 ){
   // Extraction : 1.3     // Amplification : 28.0     // Induction : 1.0  
    run_list = {747, 779, 783};
  }
  else if( scan_type == "Drift" and scan_num == 2 ){
   // Extraction : 1.2     // Amplification : 28.0     // Induction : 1.0  
    run_list = {774, 775, 776, 778};
  }
  else if( scan_type == "Drift" and scan_num == 3 ){
   // Extraction : 1.4     // Amplification : 28.0     // Induction : 1.0  
    run_list = {780, 784, 781};
  }
  else if( scan_type == "Drift" and scan_num == 4 ){
   // Extraction : 1.5     // Amplification : 28.0     // Induction : 1.0  
    run_list = {782, 785};
  }
  else if( scan_type == "Drift" and scan_num == 5 ){
   // Extraction : 1.6     // Amplification : 28.0     // Induction : 1.0  
    run_list = {786, 797};
  }
  else if( scan_type == "Drift" and scan_num == 6 ){
   // Extraction : 1.7     // Amplification : 28.0     // Induction : 1.0  
    run_list = {787, 798};
  }
  else if( scan_type == "Drift" and scan_num == 7 ){
   // Extraction : 1.8     // Amplification : 28.0     // Induction : 1.0  
    run_list = {788, 799};
  }
  else if( scan_type == "Drift" and scan_num == 8 ){
   // Extraction : 1.9     // Amplification : 28.0     // Induction : 1.5  
    run_list = {836, 840, 841, 842, 843};
  }
  else if( scan_type == "Drift" and scan_num == 9 ){
   // Extraction : 2.0     // Amplification : 28.0     // Induction : 1.5  
    run_list = {837, 838};
  }
  else if( scan_type == "Drift" and scan_num == 10 ){
   // Extraction : 1.0     // Amplification : 30.0     // Induction : 2.0  
    run_list = {1167, 1193, 1197};
  }
  else if( scan_type == "Drift" and scan_num == 11 ){
   // Extraction : 1.6     // Amplification : 23.0     // Induction : 2.5  
    run_list = {1181, 1189};
  }

  //Generated with python macro (all runs, 50% drift tolerence)
  // Extraction :
  else if( scan_type == "Extraction" and scan_num == 1 ){
   // Drift : 0.321     // Amplification : 26.0     // Induction : 1.0  
    run_list = {743, 744};
  }
  else if( scan_type == "Extraction" and scan_num == 2 ){
   // Drift : 0.347     // Amplification : 28.0     // Induction : 1.0  
   //runs 787 and 798 have same extr field but quite different gain-> One did not have min_dQ, can change fit
    run_list = {747, 769, 770, 774, 776, 778, 780, 781, 782, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 779, 783};
  }
  else if( scan_type == "Extraction" and scan_num == 3 ){
   // Drift : 0.31     // Amplification : 29.0     // Induction : 1.0  
    run_list = {766, 767};
  }
  else if( scan_type == "Extraction" and scan_num == 4 ){
   // Drift : 0.523     // Amplification : 28.0     // Induction : 1.5  
    run_list = {833, 834, 835, 836, 837, 838, 840, 841, 842, 843};
  }
  else if( scan_type == "Extraction" and scan_num == 5 ){
   // Drift : 0.509     // Amplification : 27.0     // Induction : 1.0  
    run_list = {1000, 1002, 1003, 1004};
  }
//  else if( scan_type == "Extraction" and scan_num == 6 ){
//   // Drift : 0.207     // Amplification : 30.0     // Induction : 2.5  
//    run_list = {1165, 1166};
//  }
  else if( scan_type == "Extraction" and scan_num == 7 ){
   // Drift : 0.257     // Amplification : 23.0     // Induction : 5.0  
    run_list = {1172, 1173};
  }
  else if( scan_type == "Extraction" and scan_num == 8 ){
   // Drift : 0.36     // Amplification : 23.0     // Induction : 2.5  
    run_list = {1178, 1180, 1181, 1189};
  }

  // Amplification :
  else if( scan_type == "Amplification" and scan_num == 1 ){
   // Drift : 0.313     // Extraction : 1.0     // Induction : 1.0  
    run_list = {741, 769, 770};
  }
  else if( scan_type == "Amplification" and scan_num == 2 ){
   // Drift : 0.33     // Extraction : 1.2     // Induction : 1.0  
    run_list = {744, 745, 774, 776, 778};
  }
  else if( scan_type == "Amplification" and scan_num == 3 ){
   // Drift : 0.31     // Extraction : 0.9     // Induction : 1.0  
    run_list = {766, 768};
  }
  else if( scan_type == "Amplification" and scan_num == 4 ){
   // Drift : 0.509     // Extraction : 2.1     // Induction : 1.0  
    run_list = {791, 1005, 1006, 802};
  }
  else if( scan_type == "Amplification" and scan_num == 5 ){
   // Drift : 0.508     // Extraction : 2.2     // Induction : 1.0  
    run_list = {792, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 993, 994, 995, 996, 997, 998, 999, 1000, 803};
  }
  else if( scan_type == "Amplification" and scan_num == 6 ){
   // Drift : 0.36     // Extraction : 1.2     // Induction : 2.5  
    run_list = {1178, 1191};
  }
  else if( scan_type == "Amplification" and scan_num == 7 ){
   // Drift : 0.358     // Extraction : 1.4     // Induction : 2.5  
    run_list = {1180, 1190};
  }
  //from Elog
  else if(scan_type == "Amplification" and scan_num == 100){
    run_list = {/*1165, 1166, 1167, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1180, 1181, */1182, 1183, 1187, 1188, 1189, 1190, 1191, /*1192, 1193,*/ 1194, 1195, 1196, 1197, /*1198, 1199,*/ 840};
  }

  // Induction :
  else if( scan_type == "Induction" and scan_num == 1 ){
   // Drift : 0.454     // Extraction : 1.6     // Amplification : 28.0  
    run_list = {786, 833, 797};
  }
  else if( scan_type == "Induction" and scan_num == 2 ){
   // Drift : 0.453     // Extraction : 1.7     // Amplification : 28.0  
    run_list = {787, 834, 798};
  }
  else if( scan_type == "Induction" and scan_num == 3 ){
   // Drift : 0.452     // Extraction : 1.8     // Amplification : 28.0  
    run_list = {788, 835, 799};
  }
  else if( scan_type == "Induction" and scan_num == 4 ){
   // Drift : 0.511     // Extraction : 1.9     // Amplification : 28.0  
    run_list = {789, 836, 840, 841, 842, 843, 800};
  }
  else if( scan_type == "Induction" and scan_num == 5 ){
   // Drift : 0.51     // Extraction : 2.0     // Amplification : 28.0  
    run_list = {790, 837, 838, 801};
  }
  
  //All
  else if( scan_type == "All" and scan_num == 1 ){
    run_list = {741, 769, 770};
  }
  else if( scan_type == "All" and scan_num == 2 ){
    run_list = {743};
  }
  else if( scan_type == "All" and scan_num == 3 ){
    run_list = {744, 745, 774, 776, 778};
  }
  else if( scan_type == "All" and scan_num == 4 ){
    run_list = {747, 779};
  }
  else if( scan_type == "All" and scan_num == 5 ){
    run_list = {765};
  }
  else if( scan_type == "All" and scan_num == 6 ){
    run_list = {766, 768};
  }
  else if( scan_type == "All" and scan_num == 7 ){
    run_list = {767};
  }
  else if( scan_type == "All" and scan_num == 8 ){
    run_list = {775};
  }
  else if( scan_type == "All" and scan_num == 9 ){
    run_list = {780, 781};
  }
  else if( scan_type == "All" and scan_num == 10 ){
    run_list = {782};
  }
  else if( scan_type == "All" and scan_num == 11 ){
    run_list = {783};
  }
  else if( scan_type == "All" and scan_num == 12 ){
    run_list = {784};
  }
  else if( scan_type == "All" and scan_num == 13 ){
    run_list = {785};
  }
  else if( scan_type == "All" and scan_num == 14 ){
    run_list = {786};
  }
  else if( scan_type == "All" and scan_num == 15 ){
    run_list = {787};
  }
  else if( scan_type == "All" and scan_num == 16 ){
    run_list = {788};
  }
  else if( scan_type == "All" and scan_num == 17 ){
    run_list = {789, 800};
  }
  else if( scan_type == "All" and scan_num == 18 ){
    run_list = {790, 801};
  }
  else if( scan_type == "All" and scan_num == 19 ){
    run_list = {791, 802, 1005, 1006};
  }
  else if( scan_type == "All" and scan_num == 20 ){
    run_list = {792, 803, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 993, 994, 995, 996, 997, 998, 999, 1000};
  }
  else if( scan_type == "All" and scan_num == 21 ){
    run_list = {793, 804, 805, 806, 807};
  }
  else if( scan_type == "All" and scan_num == 22 ){
    run_list = {794};
  }
  else if( scan_type == "All" and scan_num == 23 ){
    run_list = {795};
  }
  else if( scan_type == "All" and scan_num == 24 ){
    run_list = {796};
  }
  else if( scan_type == "All" and scan_num == 25 ){
    run_list = {797};
  }
  else if( scan_type == "All" and scan_num == 26 ){
    run_list = {798};
  }
  else if( scan_type == "All" and scan_num == 27 ){
    run_list = {799};
  }
  else if( scan_type == "All" and scan_num == 28 ){
    run_list = {833};
  }
  else if( scan_type == "All" and scan_num == 29 ){
    run_list = {834};
  }
  else if( scan_type == "All" and scan_num == 30 ){
    run_list = {835};
  }
  else if( scan_type == "All" and scan_num == 31 ){
    run_list = {836, 840, 841, 842, 843};
  }
  else if( scan_type == "All" and scan_num == 32 ){
    run_list = {837};
  }
  else if( scan_type == "All" and scan_num == 33 ){
    run_list = {838};
  }
  else if( scan_type == "All" and scan_num == 34 ){
    run_list = {1002, 1003, 1004};
  }
  else if( scan_type == "All" and scan_num == 35 ){
    run_list = {1007, 1008, 1009, 1010, 1011};
  }
  else if( scan_type == "All" and scan_num == 36 ){
    run_list = {1012, 1013, 1014, 1016, 1035, 1036, 1037, 1038};
  }
  else if( scan_type == "All" and scan_num == 37 ){
    run_list = {1039, 1040};
  }
  else if( scan_type == "All" and scan_num == 38 ){
    run_list = {1165};
  }
//  else if( scan_type == "All" and scan_num == 39 ){
//    run_list = {1166};
//  }
  else if( scan_type == "All" and scan_num == 40 ){
    run_list = {1167};
  }
  else if( scan_type == "All" and scan_num == 41 ){
    run_list = {1172};
  }
  else if( scan_type == "All" and scan_num == 42 ){
    run_list = {1173};
  }
  else if( scan_type == "All" and scan_num == 43 ){
    run_list = {1174};
  }
  else if( scan_type == "All" and scan_num == 44 ){
    run_list = {1175};
  }
  else if( scan_type == "All" and scan_num == 45 ){
    run_list = {1176};
  }
  else if( scan_type == "All" and scan_num == 46 ){
    run_list = {1177};
  }
  else if( scan_type == "All" and scan_num == 47 ){
    run_list = {1178};
  }
  else if( scan_type == "All" and scan_num == 48 ){
    run_list = {1180};
  }
  else if( scan_type == "All" and scan_num == 49 ){
    run_list = {1181};
  }
  else if( scan_type == "All" and scan_num == 50 ){
    run_list = {1182};
  }
  else if( scan_type == "All" and scan_num == 51 ){
    run_list = {1183};
  }
  else if( scan_type == "All" and scan_num == 52 ){
    run_list = {1187};
  }
  else if( scan_type == "All" and scan_num == 53 ){
    run_list = {1188};
  }
  else if( scan_type == "All" and scan_num == 54 ){
    run_list = {1189};
  }
  else if( scan_type == "All" and scan_num == 55 ){
    run_list = {1190};
  }
  else if( scan_type == "All" and scan_num == 56 ){
    run_list = {1191};
  }
  else if( scan_type == "All" and scan_num == 57 ){
    run_list = {1192};
  }
  else if( scan_type == "All" and scan_num == 58 ){
    run_list = {1193, 1197};
  }
  else if( scan_type == "All" and scan_num == 59 ){
    run_list = {1194};
  }
  else if( scan_type == "All" and scan_num == 60 ){
    run_list = {1195};
  }
  else if( scan_type == "All" and scan_num == 61 ){
    run_list = {1196};
  }
  else if( scan_type == "All" and scan_num == 62 ){
    run_list = {1198};
  }
  else if( scan_type == "All" and scan_num == 63 ){
    run_list = {1199};
  }
  
  //Identical runs:
  else if( scan_type == "Identical" and scan_num == 1 ){
   // Drift : 0.347     // Extraction : 1.3     // Amplification : 28.0     // Induction : 1.0  
    run_list = {747, 779, 783};
  }
  else if( scan_type == "Identical" and scan_num == 2 ){
   // Drift : 0.31     // Extraction : 1.0     // Amplification : 28.0     // Induction : 1.0  
    run_list = {769, 770};
  }
  else if( scan_type == "Identical" and scan_num == 3 ){
   // Drift : 0.358     // Extraction : 1.2     // Amplification : 28.0     // Induction : 1.0  
    run_list = {774, 776, 778};
  }
  else if( scan_type == "Identical" and scan_num == 4 ){
   // Drift : 0.356     // Extraction : 1.4     // Amplification : 28.0     // Induction : 1.0  
    run_list = {780, 781, 784};
  }
  else if( scan_type == "Identical" and scan_num == 5 ){
   // Drift : 0.355     // Extraction : 1.5     // Amplification : 28.0     // Induction : 1.0  
    run_list = {782, 785};
  }
  else if( scan_type == "Identical" and scan_num == 6 ){
   // Drift : 0.454     // Extraction : 1.6     // Amplification : 28.0     // Induction : 1.0  
    run_list = {786, 797};
  }
  else if( scan_type == "Identical" and scan_num == 7 ){
   // Drift : 0.453     // Extraction : 1.7     // Amplification : 28.0     // Induction : 1.0  
    run_list = {787, 798};
  }
  else if( scan_type == "Identical" and scan_num == 8 ){
   // Drift : 0.452     // Extraction : 1.8     // Amplification : 28.0     // Induction : 1.0  
    run_list = {788, 799};
  }
  else if( scan_type == "Identical" and scan_num == 9 ){
   // Drift : 0.511     // Extraction : 1.9     // Amplification : 28.0     // Induction : 1.0  
    run_list = {789, 800};
  }
  else if( scan_type == "Identical" and scan_num == 10 ){
   // Drift : 0.51     // Extraction : 2.0     // Amplification : 28.0     // Induction : 1.0  
    run_list = {790, 801};
  }
  else if( scan_type == "Identical" and scan_num == 11 ){
   // Drift : 0.509     // Extraction : 2.1     // Amplification : 28.0     // Induction : 1.0  
    run_list = {791, 802};
  }
  else if( scan_type == "Identical" and scan_num == 12 ){
   // Drift : 0.508     // Extraction : 2.2     // Amplification : 28.0     // Induction : 1.0  
    run_list = {792, 803};
  }
  else if( scan_type == "Identical" and scan_num == 13 ){
   // Drift : 0.507     // Extraction : 2.3     // Amplification : 28.0     // Induction : 1.0  
    run_list = {793, 804, 805, 806, 807};
  }
  else if( scan_type == "Identical" and scan_num == 14 ){
   // Drift : 0.52     // Extraction : 1.9     // Amplification : 28.0     // Induction : 1.5  
    run_list = {836, 840, 841, 842, 843};
  }
  else if( scan_type == "Identical" and scan_num == 15 ){
   // Drift : 0.519     // Extraction : 2.0     // Amplification : 28.0     // Induction : 1.5  
    run_list = {837, 838};
  }
  else if( scan_type == "Identical" and scan_num == 16 ){
   // Drift : 0.511     // Extraction : 2.2     // Amplification : 25.0     // Induction : 1.0  
    run_list = {981, 982, 983, 984, 985, 986, 987, 988, 989};
  }
  else if( scan_type == "Identical" and scan_num == 17 ){
   // Drift : 0.5105     // Extraction : 2.2     // Amplification : 25.5     // Induction : 1.0  
    run_list = {990, 991, 993, 994};
  }
  else if( scan_type == "Identical" and scan_num == 18 ){
   // Drift : 0.51     // Extraction : 2.2     // Amplification : 26.0     // Induction : 1.0  
    run_list = {995, 996};
  }
  else if( scan_type == "Identical" and scan_num == 19 ){
   // Drift : 0.5095     // Extraction : 2.2     // Amplification : 26.5     // Induction : 1.0  
    run_list = {997, 998, 999};
  }
  else if( scan_type == "Identical" and scan_num == 20 ){
   // Drift : 0.5095     // Extraction : 2.15     // Amplification : 27.0     // Induction : 1.0  
    run_list = {1002, 1003, 1004};
  }
  else if( scan_type == "Identical" and scan_num == 21 ){
   // Drift : 0.5095     // Extraction : 2.1     // Amplification : 27.5     // Induction : 1.0  
    run_list = {1005, 1006};
  }
  else if( scan_type == "Identical" and scan_num == 22 ){
   // Drift : 0.51     // Extraction : 2.05     // Amplification : 27.0     // Induction : 1.25  
    run_list = {1007, 1008, 1009, 1010, 1011};
  }
  else if( scan_type == "Identical" and scan_num == 23 ){
   // Drift : 0.51     // Extraction : 2.0     // Amplification : 27.5     // Induction : 1.25  
    run_list = {1012, 1013, 1014, 1016, 1035, 1036, 1037, 1038};
  }
  else if( scan_type == "Identical" and scan_num == 24 ){
   // Drift : 0.51     // Extraction : 1.95     // Amplification : 28.0     // Induction : 1.25  
    run_list = {1039, 1040};
  }
  else if( scan_type == "Identical" and scan_num == 25 ){
   // Drift : 0.356     // Extraction : 1.6     // Amplification : 23.0     // Induction : 2.5  
    run_list = {1181, 1189};
  }
  else if( scan_type == "Identical" and scan_num == 26 ){
   // Drift : 0.516     // Extraction : 1.0     // Amplification : 30.0     // Induction : 2.0  
    run_list = {1193, 1197};
  }

  else{
    #if verbose
    cout<< "Error: unknown combination of scan type and scan number (" << scan_type << ", " << scan_num << ")." << endl;
    #endif
    return false;
  }
  if(scan_type == "Identical"){
    return true;
  }
  for(auto r : run_list){
    field.push_back(runs_and_fields[r][field_to_read]);
  }
  for(auto f : runs_and_fields[run_list[0]]){
    if(f.first == field_to_read){continue;}
    const_fields_names.push_back(f.first);
    const_fields_values.push_back(f.second);
  }
  return true;
}


bool load_run_lists(){
  #if verbose
  cout << "Loading runs from " << runs_file << "..." << endl;
  #endif
  vector<string> field_names;
  fstream all_runs;
  all_runs.open(runs_file.data(), fstream::in);
  if(!all_runs.is_open()){
    cout << "ERROR : can not open all_runs.csv" << endl;
    return false;
  }
  int i = 0;
  while(!all_runs.eof()){
    char buffer[120];
    all_runs.getline(buffer, 120);
    char* tokBuffer;
    if( string(buffer).find(";") == string::npos ){break;}
    tokBuffer = strtok(buffer, ";");
    int run_tmp = 0;
    if(i == 0){
      string dummy = string(tokBuffer);
      tokBuffer = strtok(NULL, ";");
      field_names.push_back(string(tokBuffer));
      tokBuffer = strtok(NULL, ";");
      field_names.push_back(string(tokBuffer));
      tokBuffer = strtok(NULL, ";");
      field_names.push_back(string(tokBuffer));
      tokBuffer = strtok(NULL, ";");
      field_names.push_back(string(tokBuffer));
      tokBuffer = strtok(NULL, ";");
      field_names.push_back(string(tokBuffer));
      i++; 
      continue;
    }
    run_tmp = atoi(tokBuffer);
    tokBuffer = strtok(NULL, ";");
    runs_and_fields[run_tmp][field_names[0]] = atof(tokBuffer)/1000.;
    tokBuffer = strtok(NULL, ";");
    runs_and_fields[run_tmp][field_names[1]] = atof(tokBuffer)/1000.;
    tokBuffer = strtok(NULL, ";");
    runs_and_fields[run_tmp][field_names[2]] = atof(tokBuffer)/1000.;
    tokBuffer = strtok(NULL, ";");
    runs_and_fields[run_tmp][field_names[3]] = atof(tokBuffer)/1000.;
    tokBuffer = strtok(NULL, ";");
    runs_triggers[run_tmp] = tokBuffer;
  }
  all_runs.close();
  if(runs_and_fields.size() == 0){
    cout << "ERROR loading runs : empty run list" << endl;
    return false;
  }
  cout << "...done" << endl;
  return true;
}

bool load_extr_eff_simu_graphs(){
  string path_ExtrEff_vs_LemExtr = path_311analysis + "ExtrEff_vs_LemExtr.root";
  if(!ExistTest(path_ExtrEff_vs_LemExtr)){
    cout << "ERROR: could not find " << path_ExtrEff_vs_LemExtr << endl;
    return false;
  }
  TFile ifile_ExtrEff_vs_LemExtr(path_ExtrEff_vs_LemExtr.data(), "READ");
  
  TCanvas *dummycan0;
  ifile_ExtrEff_vs_LemExtr.GetObject("can_eff_vs_LemExtr", dummycan0);
  TList *primitives0 = dummycan0->GetListOfPrimitives();
  TGraph2D *graph0 = (TGraph2D*)primitives0->At(0);
  h_ExtrEff_vs_LemExtr = *graph0->GetHistogram();
  
  if(h_ExtrEff_vs_LemExtr.GetEntries() == 0){
    cout << "Error while loading graph ExtrEff_vs_LemExtr" << endl;
    return false;
  }
  
  return true;
}
  
bool load_ind_eff_simu_graphs(){
  string path_IndEff_vs_LemInd = path_311analysis + "IndEff_vs_LemInd.root";
  if(!ExistTest(path_IndEff_vs_LemInd)){
    cout << "ERROR: could not find " << path_IndEff_vs_LemInd << endl;
    return false;
  }
  TFile ifile_IndEff_vs_LemInd(path_IndEff_vs_LemInd.data(), "READ");
  
  TCanvas *dummycan1;
  ifile_IndEff_vs_LemInd.GetObject("can_eff_vs_LemInd", dummycan1);
  TList *primitives1 = dummycan1->GetListOfPrimitives();
  TGraph2D *graph1 = (TGraph2D*)primitives1->At(0);
  h_IndEff_vs_LemInd = *graph1->GetHistogram();

  if(h_IndEff_vs_LemInd.GetEntries() == 0){
    cout << "Error while loading graph IndEff_vs_LemInd" << endl;
    return false;
  }
  return true;
}

TMyFileHeader load_run_header(int run, bool clean){
  string filepath = runs_headers + to_string(run) + ".root";
  TMyFileHeader header = TMyFileHeader();
  if(!ExistTest(filepath)){if(!save_runs_headers({run})){return header;}}
  TFile fileheader(filepath.data(),"READ");
  if(!fileheader.IsOpen()){
    cout << "ERROR: can not find header file for run " << run << endl;
    return header;
  }
  header = *(TMyFileHeader*)((TKey*)fileheader.GetListOfKeys()->At(0))->ReadObj();
  fileheader.Close();
  if(clean){
    header.Clean();
  }
  return header;
}

void load_fit_3L(){
  vector<double> fields = {29,30,31,32,33,34,34.5,35};
  vector<double> fields_err = {};
  vector<double> gains = {7.8,11.7,18.5,30.4,51.4,89.6,122.1,166.3};
  vector<double> gains_errors = {};
  for(auto g : gains){
    gains_errors.push_back(g*0.02);
    fields_err.push_back(0);
  }
  Gain_graph_3L = TGraphErrors(fields.size(), &fields[0], &gains[0], &fields_err[0], &gains_errors[0]);
  TF1 *func = new TF1("Theoretical_Gain","[0]*exp([1]*0.1*exp(-[2]/x))",0,40);
  func->SetParameter(0,1);
  func->SetParName(0,"T");
  func->SetParLimits(0,0,1);
  func->SetParameter(1,4000);
  func->SetParName(1,"A*rho");
  func->SetParameter(2,200);
  func->SetParName(2,"B*rho");
  Gain_graph_3L.Fit(func);
//  Gain_graph.Draw("A*");
//  gPad->SaveAs(string(path_311analysis+"gain/3L.png").data());
  
  Arho_from_3L = make_pair(func->GetParameter(1),func->GetParError(1)); //pressure 980mbar, temperature 87K
  Brho_from_3L = make_pair(func->GetParameter(2),func->GetParError(2));

  return;
}


void load_Gushin_Eff(){
  for(int i = 0; i < GushinEff.GetN(); i++){
    double field, eff;
    GushinEff.GetPoint(i,field,eff);
    h_GushinEff.SetBinContent(i+1, eff);
  }
  return;
}

void load_gain_eff_corrections(){
  #if verbose
  cout << "Initialising hard-coded efficiency corrections..." << endl;
  #endif
  if(h_ExtrEff_vs_LemExtr.GetEntries() == 0){
    if(!load_extr_eff_simu_graphs()){return;}
  }
  if(h_IndEff_vs_LemInd.GetEntries() == 0){
    if(!load_ind_eff_simu_graphs()){return;}
  }
  if(h_GushinEff.GetEntries() == 0){load_Gushin_Eff();}
  gain_corrections[1178] = get_eff(1.2)*.5;
  gain_corrections[1180] = get_eff(1.4)*.5;
  gain_corrections[1181] = get_eff(1.6)*.5;
  gain_corrections[1182] = get_eff(1.8)*.5;
  gain_corrections[1183] = get_eff(2.)*.5;
  gain_corrections[1187] = get_eff(2.2)*.45;
  gain_corrections[1188] = get_eff(2.4)*.4;
  gain_corrections[1189] = get_eff(1.6)*.5;
  gain_corrections[1190] = get_eff(1.4)*.5;
  gain_corrections[1191] = get_eff(1.2)*.55;
  gain_corrections[1192] = get_eff(1.)*.6;
  gain_corrections[1193] = get_eff(1.)*.6;
  gain_corrections[1194] = get_eff(1.)*.6;
  gain_corrections[1195] = get_eff(1.1)*.6;
  gain_corrections[1196] = get_eff(1.1)*.6;
  gain_corrections[1197] = get_eff(1.)*.6;
  gain_corrections[1198] = get_eff(1.1)*.6;
  gain_corrections[1199] = get_eff(1.1)*.6;
  return;
}

bool load_highway(int run){

  #if verbose
  cout << "Loading highway tracks selection output for run " << run << "..." << endl;
  #endif
  string path = highway_output + "Run" + to_string(run) + "/";
  if(!ExistTest(path)){
    tracks_selected_by_highway.push_back(-1.);
    return true;
  }
  string wildcard_path = path+"*";
  vector<string> ifiles;
  for( auto file : glob(wildcard_path) ){
    ifiles.push_back(path+file);
  }
  if(ifiles.size() == 0){
    tracks_selected_by_highway.push_back(-1.);
    return true;
  }
  
  
  for (auto file : ifiles){
    fstream ifile;
    ifile.open(file.data(), fstream::in);
    if(!ifile.is_open()){
      cout << "ERROR : can not open file " << file.data() << endl;
      return false;
    }
    int i = 0;
    
    while(!ifile.eof()){
      char buffer[1000];
      ifile.getline(buffer, 1000);
      if(i == 0){i++; continue;}
      char* tokBuffer;
      if( string(buffer).find(";") == string::npos ){break;}
      tokBuffer = strtok(buffer, ";"); tokBuffer = strtok(NULL, ";"); tokBuffer = strtok(NULL, ";");
      long track_ID = 1000*(1+atoi(tokBuffer));
      tokBuffer = strtok(NULL, ";");
      track_ID += atoi(tokBuffer);
      tracks_selected_by_highway.push_back(track_ID);
    }
    ifile.close();
  }
  
  #if verbose
  cout << "...Done. Found " << tracks_selected_by_highway.size() << " tracks." << endl;
  #endif
  return true;
}


void load_force_mpv(){
  force_mpv[1182] = make_pair(2.220945, 2.053639);
  force_mpv[1183] = make_pair(2.184064, 2.062411);
  force_mpv[1187] = make_pair(1.921035, 2.049250);
  force_mpv[1188] = make_pair(2.033270, 2.055467);
  force_mpv[1189] = make_pair(2.876688, 2.564887);
  force_mpv[1190] = make_pair(2.997272, 2.662352);
  force_mpv[1191] = make_pair(3.007459, 2.811454);
  force_mpv[1194] = make_pair(3.574564, 3.406592);
  force_mpv[1195] = make_pair(2.421819, 2.656836);
  force_mpv[1196] = make_pair(3.564599, 2.998725);
  force_mpv[1197] = make_pair(3.641836, 3.092138);
  return;
}

bool load_gr_and_h_rho(int run){
  string ifilepath = slow_control + to_string(run) + "/pressure_and_temp/temperature_cryostat.root";
  if(!ExistTest(ifilepath)){
    #if verbose
    cout << "Plotting pressure, temp and density graphs for run " << run <<"..." << endl;
    #endif
    if(!pressure(to_string(run))){return false;}
  }
  TFile ifile(ifilepath.data(), "READ");
  if(!ifile.IsOpen()){
    cout << "ERROR in load_gr_and_h_rho: can not open file " << string(slow_control + to_string(run) + "/pressure_and_temp/temperature_cryostat.root").data() << endl;
    return false;
  }
  TGraph *dummygraph = 0;
  TH1D *dummyhisto = 0;
  ifile.GetObject("graph_PE0006OverTE0056",dummygraph);
  ifile.GetObject("histo_PE0006OverTE0056",dummyhisto);
  if(dummygraph == 0){
    cout << "ERROR: can not get graph_PE0006OverTE0056 in file " << ifile.GetName() << endl;
    return false;
  }
  Gr_Rho = TGraph(*dummygraph);
  H_Rho = TH1D(*dummyhisto);
  return true;
}

int find_lem(double y, double z){
  //for the 311 geometry retun the lem nubmer associated to a certain geometry
  int lem =-99999;
  int lem_in_module = 4; //num of minimal modules of 4 lems
  if(y == tpc_boundaries[3]){y = y-1e-5;}
  if(z == tpc_boundaries[5]){z = z-1e-5;}
 //fist of all check if the set of coodinate is valid (tolerance on boundaries)
  if( fabs(z) >= tpc_boundaries[5] or fabs(y) >= tpc_boundaries[3]){
//    #if verbose
//    cout << "z " << z << " and y " << y <<" are unknown coordinates" << endl;
//    #endif
    return lem;
  }
  
       if(z >=    0          and z < lem_size   and y >= 0         and y < lem_size ){lem =  1; }
  else if(z >=    1*lem_size and z < 2*lem_size and y >= 0         and y < lem_size ){lem =  2; }
  else if(z >=    0          and z < lem_size   and y >= -lem_size and y < 0        ){lem =  3; }
  else if(z >=    1*lem_size and z < 2*lem_size and y >= -lem_size and y < 0        ){lem =  4; }
  else if(z >=    2*lem_size and z < 3*lem_size and y >= 0         and y < lem_size ){lem =  5; }
  else if(z >=    3*lem_size and z < 4*lem_size and y >= 0         and y < lem_size ){lem =  6; }
  else if(z >=    2*lem_size and z < 3*lem_size and y >= -lem_size and y < 0        ){lem =  7; }
  else if(z >=    3*lem_size and z < 4*lem_size and y >= -lem_size and y < 0        ){lem =  8; }
  else if(z >=    4*lem_size and z < 5*lem_size and y >= 0         and y < lem_size ){lem =  9; }
  else if(z >=    5*lem_size and z < 6*lem_size and y >= 0         and y < lem_size ){lem = 10; }
  else if(z >=    4*lem_size and z < 5*lem_size and y >= -lem_size and y < 0        ){lem = 11; }
  else if(z >=    5*lem_size and z < 6*lem_size and y >= -lem_size and y < 0        ){lem = 12; }

  return lem;
}


pair<int,int> find_channels(double y, double z){
  if(y == tpc_boundaries[3]){y = y-1e-5;}
  if(z == tpc_boundaries[5]){z = z-1e-5;}
  int ch0 = (int)((y-tpc_boundaries[2])/pitch);
  int ch1 = N_Ch_0 + (int)((z-tpc_boundaries[4])/pitch);
  return make_pair(ch0,ch1);
}
bool IsGoodChan(int ch0, int ch1){
  if( find(bad_channels.begin(), bad_channels.end(), ch0) != bad_channels.end() or find(bad_channels.begin(), bad_channels.end(), ch1) != bad_channels.end() ){return false;}
  else if( ch0 < 0 or ch0 > N_Ch_0 - 1 ){return false;}
  else if( ch1 < N_Ch_0 or ch1 > N_Ch_1 +N_Ch_0 - 1 ){return false;}
  return true;
}
bool IsGoodChan(double y, double z){
  pair<int,int> chans = find_channels(y, z);
  return IsGoodChan(chans.first, chans.second);
}
bool IsGoodChan(hit myhit){
  return IsGoodChan(myhit.sp_y, myhit.sp_z);
}

vector<int> find_YZ(int lem){
  int y =-100;
  int z = -100;
  vector<int> to_return = {y,z};
  
       if( lem == 1 ) { to_return[0] = 0;          to_return[1] = 0;          }
  else if( lem == 2 ) { to_return[0] = 0;          to_return[1] = lem_size;   }
  else if( lem == 3 ) { to_return[0] = -lem_size;  to_return[1] = 0;          }
  else if( lem == 4 ) { to_return[0] = -lem_size;  to_return[1] = lem_size;   }
  else if( lem == 5 ) { to_return[0] = 0;          to_return[1] = 2*lem_size; }
  else if( lem == 6 ) { to_return[0] = 0;          to_return[1] = 3*lem_size; }
  else if( lem == 7 ) { to_return[0] = -lem_size;  to_return[1] = 2*lem_size; }
  else if( lem == 8 ) { to_return[0] = -lem_size;  to_return[1] = 3*lem_size; }
  else if( lem == 9 ) { to_return[0] = 0;          to_return[1] = 4*lem_size; }
  else if( lem == 10 ) { to_return[0] = 0;         to_return[1] = 5*lem_size; }
  else if( lem == 11 ) { to_return[0] = -lem_size; to_return[1] = 4*lem_size; }
  else if( lem == 12 ) { to_return[0] = -lem_size; to_return[1] = 5*lem_size; }

  return to_return;
}

bool isGood_lem( int lem ){
  //select only good lems among a list

   if (find(lems.begin(),lems.end(), lem) != lems.end())
     return true;
   else
      return false;
}

bool IsGood_run(string cut_type_and_methods, int run){
  if( bad_runs.find(cut_type_and_methods) != bad_runs.end() ){
    if( find(bad_runs[cut_type_and_methods].begin(), bad_runs[cut_type_and_methods].end(), run) != bad_runs[cut_type_and_methods].end() ){
      return false;
    }
  }
  return true;
}
bool IsGood_lem_gain(string cut_type_and_methods, int run, int lem){
  if(!IsGood_run(cut_type_and_methods, run)){return false;}
  if( bad_runs_lems.find(cut_type_and_methods) != bad_runs_lems.end() ){
    if( bad_runs_lems[cut_type_and_methods].find(run) != bad_runs_lems[cut_type_and_methods].end() ){
      if( find(bad_runs_lems[cut_type_and_methods][run].begin(), bad_runs_lems[cut_type_and_methods][run].end(), lem) != bad_runs_lems[cut_type_and_methods][run].end() ){
        return false;
      }
    }
  }
  return true;
}
bool IsGood_yz(string cut_type_and_methods, int run, pair<int,int> YZ){
  if(!IsGood_run(cut_type_and_methods, run)){return false;}
    if( bad_runs_yz.find(cut_type_and_methods) != bad_runs_yz.end() ){
      if( bad_runs_yz[cut_type_and_methods].find(run) != bad_runs_yz[cut_type_and_methods].end() ){
        if( find(bad_runs_yz[cut_type_and_methods][run].begin(), bad_runs_yz[cut_type_and_methods][run].end(), YZ) != bad_runs_yz[cut_type_and_methods][run].end() ){
        return false;
      }
    }
  }
  return true;
}

double find_projection(hit h){
  //return the correct variable projection on the view
  double coor = -99999;
  switch (h.view) {
    case 0:
      coor = h.sp_y;
      break;
    case 1:
      coor = h.sp_z;
      break;
    default:
      cout << " Unknown projection " << endl;
  }
  return coor;
}

double get_theta(track t){
  double theta = fabs(acos( ( t.end_x - t.start_x )/abs( t.length ) ))*180./TMath::Pi();
  if ( theta > 180 ){
    theta = 360-theta;
  }
  if ( theta < 90 ){
    theta = 180-theta;
  }
  return theta;
}

double get_phi(track t){
  if( abs(t.end_z - t.start_z)!=0 ){
    if( (t.end_z - t.start_z) > 0 ){
      return atan( (t.end_y - t.start_y)/(t.end_z - t.start_z) )*180./TMath::Pi();
    }
    else if( atan( (t.end_y - t.start_y)/(t.end_z - t.start_z) ) < 0 ){
      return atan( (t.end_y - t.start_y)/(t.end_z - t.start_z) )*180./TMath::Pi()+180.;
    }
    else{
      return atan( (t.end_y - t.start_y)/(t.end_z - t.start_z) )*180./TMath::Pi()-180.;
    }
  }
  else if( (t.end_y - t.start_y) > 0 ){
    return 90.;
  }
  else{
    return -90.;
  }
}

double get_reduced_phi(track t){
  if( abs(t.end_z - t.start_z)!=0 ){
    double phi = atan( fabs(t.end_y - t.start_y)/fabs(t.end_z - t.start_z) )*180./TMath::Pi();
    if (phi > 90){
      phi = 180-phi;
    }
    return phi;
  }
  else{
    return 90;
  }
}
  
vector<double> get_dss(track t){
  double ds_0 = pitch/fabs(TMath::Sin(get_theta(t))*TMath::Sin(get_phi(t)));
  double ds_1 = pitch/fabs(TMath::Sin(get_theta(t))*TMath::Cos(get_phi(t)));
  return {ds_0, ds_1};
}

vector<double> get_dss_from_angles(double theta, double phi){
  double ds_0 = pitch/fabs(TMath::Sin(theta*TMath::Pi()/180.)*TMath::Sin(phi*TMath::Pi()/180.));
  double ds_1 = pitch/fabs(TMath::Sin(theta*TMath::Pi()/180.)*TMath::Cos(phi*TMath::Pi()/180.));
  return {ds_0,ds_1};
}

//double get_theta(track t){
//  return acos( ( t.end_x - t.start_x )/abs( t.length ) )*180./TMath::Pi();
//}

//double get_phi(track t){
//  if(abs(t.end_z - t.start_z)!=0)
//    return atan( abs(t.end_y - t.start_y)/(t.end_z - t.start_z) )*180./TMath::Pi();
//  else
//    return 90.;
//}

void read_tree_Feb(TChain *rTree, vector<track> & tracks, int &tstart, int &tend, int to_read){
  //read the tree and store all the tree information in the respective variables
  const int NMaxHitsPerEvent=100000;
  const int NMaxClustersPerEvent=10000;
  const int NMaxTracksPerEvent=1000;
  const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
  int NEventsPerRun=335;

  //Load ROOT file and tree
  int NEntries = (int)rTree->GetEntries();

  //Define variables to store the data of the ROOT file
  //Metadata
  int tRun;
  int tSubrun;
  int tEventNumberInRun;
  int  tEventTimeSeconds;
//  int  tEventTimeNanoseconds;
//  char tIsData;

  //Hit variables
  int tNumberOfHits;
//  short tHit_TPC[NMaxHitsPerEvent];
//  short tHit_View[NMaxHitsPerEvent];
//  short tHit_Channel[NMaxHitsPerEvent];
//  float tHit_PeakTime[NMaxHitsPerEvent];
//  float tHit_ChargeSummedADC[NMaxHitsPerEvent];
//  float tHit_ChargeIntegral[NMaxHitsPerEvent];
//  float tHit_PeakHeight[NMaxHitsPerEvent];
//  float tHit_StartTime[NMaxHitsPerEvent];
//  float tHit_EndTime[NMaxHitsPerEvent];
//  float tHit_Width[NMaxHitsPerEvent];
  float tHit_GoodnessOfFit[NMaxHitsPerEvent];
//  short tHit_Multiplicity[NMaxHitsPerEvent];
//  short tHit_TrackID[NMaxHitsPerEvent];
//  short tHit_ClusterID[NMaxHitsPerEvent];

  //Cluster variables
//  short tNumberOfClusters;
//  short tClusterID[NMaxClustersPerEvent];
//  short tCluster_NumberOfHits[NMaxClustersPerEvent];
//  short tCluster_View[NMaxClustersPerEvent];
//  float tCluster_ChargeIntegral[NMaxHitsPerEvent];
//  short tCluster_StartChannel[NMaxHitsPerEvent];
//  short tCluster_StartTick[NMaxHitsPerEvent];
//  short tCluster_EndChannel[NMaxHitsPerEvent];
//  short tCluster_EndTick[NMaxHitsPerEvent];
//  float tCluster_StartCharge[NMaxHitsPerEvent];
//  float tCluster_StartAngle[NMaxHitsPerEvent];
//  float tCluster_EndCharge[NMaxHitsPerEvent];
//  float tCluster_EndAngle[NMaxHitsPerEvent];

  //Track variables
  short tNumberOfTracks;
//  short tTrackID[NMaxTracksPerEventTimesNViews];
  short tTrack_NumberOfHits[NMaxTracksPerEventTimesNViews];
  float tTrack_Length[NMaxTracksPerEventTimesNViews];

  float tTrack_StartPoint_X[NMaxTracksPerEvent];
  float tTrack_StartPoint_Y[NMaxTracksPerEvent];
  float tTrack_StartPoint_Z[NMaxTracksPerEvent];
//  float tTrack_StartPoint_DistanceToBoundary[NMaxTracksPerEvent];
  float tTrack_EndPoint_X[NMaxTracksPerEvent];
  float tTrack_EndPoint_Y[NMaxTracksPerEvent];
  float tTrack_EndPoint_Z[NMaxTracksPerEvent];
//  float tTrack_EndPoint_DistanceToBoundary[NMaxTracksPerEvent];
  
  float tTrack_StartDirection_Theta[NMaxTracksPerEvent];
  float tTrack_StartDirection_Phi[NMaxTracksPerEvent];
//  float tTrack_StartDirection_X[NMaxTracksPerEvent];
//  float tTrack_StartDirection_Y[NMaxTracksPerEvent];
//  float tTrack_StartDirection_Z[NMaxTracksPerEvent];

  float tTrack_EndDirection_Theta[NMaxTracksPerEvent];
  float tTrack_EndDirection_Phi[NMaxTracksPerEvent];
//  float tTrack_EndDirection_X[NMaxTracksPerEvent];
//  float tTrack_EndDirection_Y[NMaxTracksPerEvent];
//  float tTrack_EndDirection_Z[NMaxTracksPerEvent];

//  float tTrack_PitchInViews[NMaxTracksPerEvent];
  short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2];

  //Track hit variables
  float tTrack_Hit_X[NMaxHitsPerEvent];
  float tTrack_Hit_Y[NMaxHitsPerEvent];
  float tTrack_Hit_Z[NMaxHitsPerEvent];
  float tTrack_dx_LocalTrackDirection[NMaxHitsPerEvent];
  float tTrack_dx_3DPosition[NMaxHitsPerEvent];
//  short tTrack_Hit_TPC[NMaxHitsPerEvent];
  short tTrack_Hit_View[NMaxHitsPerEvent];
//  short tTrack_Hit_Channel[NMaxHitsPerEvent];
//  float tTrack_Hit_PeakTime[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];
//  float tTrack_Hit_PeakHeight[NMaxHitsPerEvent];
//  float tTrack_Hit_StartTime[NMaxHitsPerEvent];
//  float tTrack_Hit_EndTime[NMaxHitsPerEvent];
//  float tTrack_Hit_Width[NMaxHitsPerEvent];
  float tTrack_Hit_GoodnessOfFit[NMaxHitsPerEvent];
//  short tTrack_Hit_Multiplicity[NMaxHitsPerEvent];

  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("Run",&tRun);
  rTree->SetBranchAddress("Subrun",&tSubrun);
  rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  rTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
//  rTree->SetBranchAddress("EventTimeNanoseconds",&tEventTimeNanoseconds);
//  rTree->SetBranchAddress("IsData",&tIsData);

//  //Hit variables
  rTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);
//  rTree->SetBranchAddress("Hit_TPC",&tHit_TPC);
//  rTree->SetBranchAddress("Hit_View",&tHit_View);
//  rTree->SetBranchAddress("Hit_Channel",&tHit_Channel);
//  rTree->SetBranchAddress("Hit_ChargeSummedADC",&tHit_ChargeSummedADC);
//  rTree->SetBranchAddress("Hit_ChargeIntegral",&tHit_ChargeIntegral);
//  rTree->SetBranchAddress("Hit_PeakHeight",&tHit_PeakHeight);
//  rTree->SetBranchAddress("Hit_PeakTime",&tHit_PeakTime);
//  rTree->SetBranchAddress("Hit_StartTime",&tHit_StartTime);
//  rTree->SetBranchAddress("Hit_EndTime",&tHit_EndTime);
//  rTree->SetBranchAddress("Hit_Width",&tHit_Width);
  rTree->SetBranchAddress("Hit_GoodnessOfFit",&tHit_GoodnessOfFit);
//  rTree->SetBranchAddress("Hit_Multiplicity",&tHit_Multiplicity);
//  rTree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);
//  rTree->SetBranchAddress("Hit_ClusterID",&tHit_ClusterID);

//  //Cluster variables
//  rTree->SetBranchAddress("NumberOfClusters",&tNumberOfClusters);
//  rTree->SetBranchAddress("ClusterID",&tClusterID);
//  rTree->SetBranchAddress("Cluster_View",&tCluster_View);
//  rTree->SetBranchAddress("Cluster_NumberOfHits",&tCluster_NumberOfHits);
//  rTree->SetBranchAddress("Cluster_ChargeIntegral",&tCluster_ChargeIntegral);
//  rTree->SetBranchAddress("Cluster_StartChannel",&tCluster_StartChannel);
//  rTree->SetBranchAddress("Cluster_StartTick",&tCluster_StartTick);
//  rTree->SetBranchAddress("Cluster_EndChannel",&tCluster_EndChannel);
//  rTree->SetBranchAddress("Cluster_EndTick",&tCluster_EndTick);
//  rTree->SetBranchAddress("Cluster_StartCharge",&tCluster_StartCharge);
//  rTree->SetBranchAddress("Cluster_StartAngle",&tCluster_StartAngle);
//  rTree->SetBranchAddress("Cluster_EndCharge",&tCluster_EndCharge);
//  rTree->SetBranchAddress("Cluster_EndAngle",&tCluster_EndAngle);

//  //Track variables
  rTree->SetBranchAddress("NumberOfTracks_pmtrack",&tNumberOfTracks);
//  rTree->SetBranchAddress("TrackID_pmtrack",&tTrackID);
  rTree->SetBranchAddress("Track_NumberOfHits_pmtrack",&tTrack_NumberOfHits);
  rTree->SetBranchAddress("Track_Length_pmtrack",&tTrack_Length);

  rTree->SetBranchAddress("Track_StartPoint_X_pmtrack", &tTrack_StartPoint_X);
  rTree->SetBranchAddress("Track_StartPoint_Y_pmtrack", &tTrack_StartPoint_Y);
  rTree->SetBranchAddress("Track_StartPoint_Z_pmtrack", &tTrack_StartPoint_Z);
//  rTree->SetBranchAddress("Track_StartPoint_DistanceToBoundary_pmtrack", &tTrack_StartPoint_DistanceToBoundary);
  rTree->SetBranchAddress("Track_EndPoint_X_pmtrack", &tTrack_EndPoint_X);
  rTree->SetBranchAddress("Track_EndPoint_Y_pmtrack", &tTrack_EndPoint_Y);
  rTree->SetBranchAddress("Track_EndPoint_Z_pmtrack", &tTrack_EndPoint_Z);
//  rTree->SetBranchAddress("Track_EndPoint_DistanceToBoundary_pmtrack",&tTrack_EndPoint_DistanceToBoundary);

  rTree->SetBranchAddress("Track_StartDirection_Theta_pmtrack",&tTrack_StartDirection_Theta);
  rTree->SetBranchAddress("Track_StartDirection_Phi_pmtrack",&tTrack_StartDirection_Phi);
//  rTree->SetBranchAddress("Track_StartDirection_X_pmtrack", &tTrack_StartDirection_X);
//  rTree->SetBranchAddress("Track_StartDirection_Y_pmtrack", &tTrack_StartDirection_Y);
//  rTree->SetBranchAddress("Track_StartDirection_Z_pmtrack", &tTrack_StartDirection_Z);

  rTree->SetBranchAddress("Track_EndDirection_Theta_pmtrack",&tTrack_EndDirection_Theta);
  rTree->SetBranchAddress("Track_EndDirection_Phi_pmtrack",&tTrack_EndDirection_Phi);
//  rTree->SetBranchAddress("Track_EndDirection_X_pmtrack", &tTrack_EndDirection_X);
//  rTree->SetBranchAddress("Track_EndDirection_Y_pmtrack", &tTrack_EndDirection_Y);
//  rTree->SetBranchAddress("Track_EndDirection_Z_pmtrack", &tTrack_EndDirection_Z);

//  rTree->SetBranchAddress("Track_PitchInViews_pmtrack", &tTrack_PitchInViews);
  rTree->SetBranchAddress("Track_NumberOfHitsPerView_pmtrack",&tTrack_NumberOfHitsPerView);

//  //Track hit variables
  rTree->SetBranchAddress("Track_Hit_X_pmtrack", &tTrack_Hit_X);
  rTree->SetBranchAddress("Track_Hit_Y_pmtrack", &tTrack_Hit_Y);
  rTree->SetBranchAddress("Track_Hit_Z_pmtrack", &tTrack_Hit_Z);
  rTree->SetBranchAddress("Track_Hit_dx_LocalTrackDirection_pmtrack", &tTrack_dx_LocalTrackDirection);
  rTree->SetBranchAddress("Track_Hit_dx_3DPosition_pmtrack", &tTrack_dx_3DPosition);
//  rTree->SetBranchAddress("Track_Hit_TPC_pmtrack", &tTrack_Hit_TPC);
  rTree->SetBranchAddress("Track_Hit_View_pmtrack", &tTrack_Hit_View);
//  rTree->SetBranchAddress("Track_Hit_Channel_pmtrack", &tTrack_Hit_Channel);
//  rTree->SetBranchAddress("Track_Hit_PeakTime_pmtrack", &tTrack_Hit_PeakTime);
  rTree->SetBranchAddress("Track_Hit_ChargeSummedADC_pmtrack", &tTrack_Hit_ChargeSummedADC);
  rTree->SetBranchAddress("Track_Hit_ChargeIntegral_pmtrack", &tTrack_Hit_ChargeIntegral);
//  rTree->SetBranchAddress("Track_Hit_PeakHeight_pmtrack", &tTrack_Hit_PeakHeight);
//  rTree->SetBranchAddress("Track_Hit_StartTime_pmtrack", &tTrack_Hit_StartTime);
//  rTree->SetBranchAddress("Track_Hit_EndTime_pmtrack", &tTrack_Hit_EndTime);
//  rTree->SetBranchAddress("Track_Hit_Width_pmtrack", &tTrack_Hit_Width);
  rTree->SetBranchAddress("Track_Hit_GoodnessOfFit_pmtrack", &tTrack_Hit_GoodnessOfFit);
//  rTree->SetBranchAddress("Track_Hit_Multiplicity_pmtrack", &tTrack_Hit_Multiplicity);

  if( to_read == 0){
    to_read = NEntries;
  }
  double Efield = -1.;
  
  for(int i=0; i<to_read; i++){ //Event loop
  
    rTree->GetEntry(i);
    if(i==0){
      tstart = tEventTimeSeconds;
      TMyFileHeader header = load_run_header(tRun);
      if(header.GetRun() == -1){return;}
      Efield = header.GetAmplification();
    }
    if(i==to_read-1){
      tend = tEventTimeSeconds;
    }


    //initialize classes (not pointers for the moment)

    int a=0; //Need this counter to remember where we are in this array.
    for(int j=0; j<tNumberOfTracks; j++){ //Track loop
      
      track dummy_track;

      dummy_track.throughgoing = false;
      dummy_track.throughgoingX = false;
      dummy_track.throughgoingY = false;
      dummy_track.throughgoingZ = false;
      dummy_track.on_edge = false;
      dummy_track.highway_selected = false;
      dummy_track.blc = false;

      dummy_track.run = tRun;
      dummy_track.subrun = tSubrun;
      dummy_track.event = tEventNumberInRun;
      dummy_track.event_ntracks = tNumberOfTracks;
      dummy_track.id = j;
      dummy_track.start_x = tTrack_StartPoint_X[j];
      dummy_track.start_y = tTrack_StartPoint_Y[j];
      dummy_track.start_z = tTrack_StartPoint_Z[j];

      dummy_track.end_x = tTrack_EndPoint_X[j];
      dummy_track.end_y = tTrack_EndPoint_Y[j];
      dummy_track.end_z = tTrack_EndPoint_Z[j];
      
      
      if( max(dummy_track.end_x, dummy_track.start_x) > MaxX and min(dummy_track.end_x, dummy_track.start_x) < MinX ){
        dummy_track.throughgoingX = true;
      } //end if x
      else if( max(dummy_track.end_x, dummy_track.start_x) > MaxX and min(dummy_track.end_x, dummy_track.start_x) < MinX ){
        dummy_track.on_edge = true;
      }
      if( max(dummy_track.end_y, dummy_track.start_y) > MaxY and min(dummy_track.end_y, dummy_track.start_y) < MinY ){
        dummy_track.throughgoingY = true;
      } //end if z
      if( max(dummy_track.end_z, dummy_track.start_z) > MaxZ and min(dummy_track.end_z, dummy_track.start_z) < MinZ ){
        dummy_track.throughgoingZ = true;
      } //end if y
      if(dummy_track.throughgoingX or dummy_track.throughgoingY or dummy_track.throughgoingZ){
        dummy_track.throughgoing = true;
      }
      
      int lem0 = find_lem(0, dummy_track.start_z);
      int lem1 = find_lem(0, dummy_track.end_z);
      if(isGood_lem(lem0) and isGood_lem(lem1) ){dummy_track.blc = true;}
      if(highway){
        if(tracks_selected_by_highway.size() == 0){
          cout << "ERROR in Read_tree_Feb : please load highway output" << endl;
          return;
        }
        int track_ID = 1000*(1+dummy_track.event) + dummy_track.id;
        if( find(tracks_selected_by_highway.begin(), tracks_selected_by_highway.end(), track_ID) != tracks_selected_by_highway.end() ){
          dummy_track.highway_selected = true;
        }
      }
      dummy_track.start_theta = fabs(tTrack_StartDirection_Theta[j]);
      if( dummy_track.start_theta > 180 ){
        dummy_track.start_theta = 360 - dummy_track.start_theta;
      }
      if( dummy_track.start_theta < 90 ){
        dummy_track.start_theta = 180 - dummy_track.start_theta;
      }
      dummy_track.start_phi = tTrack_StartDirection_Phi[j];
//      dummy_track.start_vx = tTrack_StartDirection_X[j];
//      dummy_track.start_vy = tTrack_StartDirection_Y[j];
//      dummy_track.start_vz = tTrack_StartDirection_Z[j];
      
      dummy_track.end_theta = fabs(tTrack_EndDirection_Theta[j]);
      if( dummy_track.end_theta > 180 ){
        dummy_track.end_theta = 360 - dummy_track.end_theta;
      }
      if( dummy_track.end_theta < 90 ){
        dummy_track.end_theta = 180 - dummy_track.end_theta;
      }
      dummy_track.end_phi = tTrack_EndDirection_Phi[j];
//      dummy_track.end_vx = tTrack_EndDirection_X[j];
//      dummy_track.end_vy = tTrack_EndDirection_Y[j];
//      dummy_track.end_vz = tTrack_EndDirection_Z[j];

      dummy_track.length = tTrack_Length[j];

      dummy_track.theta = get_theta(dummy_track);
      dummy_track.phi = get_phi(dummy_track);
      
      
      for(int l=a; l < (a+tTrack_NumberOfHits[j]); l++){
        hit dummy_hits;

        dummy_hits.run = tRun;
        dummy_hits.subrun = tSubrun;
        dummy_hits.event = tEventNumberInRun;
        dummy_track.event_ntracks = tNumberOfTracks;
        dummy_hits.view = tTrack_Hit_View[l];
        dummy_hits.track_id = dummy_track.id;
        dummy_hits.dq_sum = tTrack_Hit_ChargeSummedADC[l]/ADC2CHARGE[tTrack_Hit_View[l]];
        dummy_hits.dq_integral = tTrack_Hit_ChargeIntegral[l]/ADC2CHARGE[tTrack_Hit_View[l]];
        dummy_hits.dx_3D = tTrack_dx_3DPosition[l];
        dummy_hits.dx_local = tTrack_dx_LocalTrackDirection[l];
        dummy_hits.sp_x = tTrack_Hit_X[l];
        dummy_hits.sp_y = tTrack_Hit_Y[l];
        dummy_hits.sp_z = tTrack_Hit_Z[l];
        dummy_hits.purity_correction = correct_dx_for_lifetime(50-tTrack_Hit_X[l]);
        dummy_hits.lem  = find_lem(dummy_hits.sp_y, dummy_hits.sp_z);
        dummy_hits.gain_density_correction_factor = gain_correction_for_rho(get_density_for_hit(tEventTimeSeconds, tRun), Efield);
        dummy_hits.GoF = tHit_GoodnessOfFit[l];
        
        dummy_track.hits_trk.push_back(dummy_hits);
      } //end hits
      tracks.push_back(dummy_track);
      a+=tTrack_NumberOfHits[j];
    }//end tracks
  }//Event loop
}

//Cuts *************************************************************************

int drays_mitigation(track & t){

  //mitigate the effects of delta rays recontructed with the track removing
  //consecutive high dqds hits
  //take track as input, return same track, but with delta rays contribution summed in their initial hit, all following n_consecutive hits deleted

  map<size_t, vector<double>> hits2view;

  unsigned int n_consecutive =15;
  double dq_cut = 6;
//  int initial_cut = 5;
  int drays = 0;
  bool Continue = false;

  unsigned int c=0;
  double sum=0;
  int view = 0;
//  int ii=0;

  vector<double> hit_list;
  double dq;
  double ds;
  double dqds;

  for(unsigned int hh=0; hh<t.hits_trk.size(); hh++){

    auto h = t.hits_trk.at(hh);
    if(method_ds == "3D"){ds = h.dx_3D;}
    else if(method_ds == "local"){ds = h.dx_local;}
    else{ds = h.dx_phil;}
    if(method_dQ == "sum"){dq = h.dq_sum;}
    else{dq = h.dq_integral;}

    //assuming ordering first all hits view 0 then all hits view 1

//    if(ii<initial_cut and h.view == view){ h.dq_sum=0; h.dq_integral=0; ++ii; }
//    else if( ii==initial_cut ){ view++; ii=0; }
    hit_list.push_back( dq );

  }

//  for(unsigned int hh=0; hh<hit_list.size(); hh++){

//    //0 pad consecutive hits above dq_cut
//    if(c<n_consecutive){
//      sum+=hit_list.at(hh);
//      ++c;
//    }
//    else if(c==n_consecutive){
//      if(sum>dq_cut*n_consecutive){
//        fill(hit_list.begin()+hh-n_consecutive, hit_list.begin()+hh, 0);
//        drays++;
//      }
//      c=0; sum=0;
//    }
//  }

  for(unsigned int hh=0; hh<hit_list.size(); hh++){

    //0 pad consecutive hits above dq_cut
    if(c==n_consecutive or Continue == true){
      if(sum>dq_cut*n_consecutive){
        Continue = false;
        fill(hit_list.begin()+hh-n_consecutive, hit_list.begin()+hh, 0);
        hit_list[hh-n_consecutive]=sum;
        c=0;
        sum=0;
        drays++;
      }
      else{
        Continue = true;
        sum = sum-hit_list.at(hh-n_consecutive);
      }
    }
    sum+=hit_list.at(hh);
    c++;
  }

  vector<hit> new_hits;
  for(unsigned int hh=0; hh<hit_list.size(); hh++){
    if( hit_list.at( hh ) > 0){
      new_hits.push_back(t.hits_trk.at(hh));
    }
  }
  t.hits_trk = new_hits;
//  for( unsigned int hh=0; hh<t.hits_trk.size(); hh++ ){
//    auto h = t.hits_trk.at(hh);
//    if(method_dQ == "sum"){t.hits_trk.at(hh).dq_sum = hit_list.at( hh );}
//    else{t.hits_trk.at(hh).dq_integral = hit_list.at( hh );}
//  }

  return drays;
}//end function


void read_tree_June(TChain *rTree, vector<track> & tracks, int &tstart, int &tend, int to_read){
  //read the tree and store all the tree information in the respective variables
  const int NMaxHitsPerEvent=100000;
  const int NMaxClustersPerEvent=10000;
  const int NMaxTracksPerEvent=1000;
  const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
  int NEventsPerRun=335;

  //Load ROOT file and tree
  int NEntries = (int)rTree->GetEntries();

  //Define variables to store the data of the ROOT file
  //Metadata
  int tRun;
  int tSubrun;
  int tEventNumberInRun;
  int  tEventTimeSeconds;
//  int  tEventTimeNanoseconds;
//  char tIsData;

  //Hit variables
  int tNumberOfHits;
//  short tHit_TPC[NMaxHitsPerEvent];
  short tHit_View[NMaxHitsPerEvent];
//  short tHit_Channel[NMaxHitsPerEvent];
//  float tHit_PeakTime[NMaxHitsPerEvent];
//  float tHit_ChargeSummedADC[NMaxHitsPerEvent];
//  float tHit_ChargeIntegral[NMaxHitsPerEvent];
//  float tHit_PeakHeight[NMaxHitsPerEvent];
//  float tHit_StartTime[NMaxHitsPerEvent];
//  float tHit_EndTime[NMaxHitsPerEvent];
//  float tHit_Width[NMaxHitsPerEvent];
  float tHit_GoodnessOfFit[NMaxHitsPerEvent];
//  short tHit_Multiplicity[NMaxHitsPerEvent];
  short tHit_TrackID[NMaxHitsPerEvent];
//  short tHit_ClusterID[NMaxHitsPerEvent];

  //Cluster variables
//  short tNumberOfClusters;
//  short tClusterID[NMaxClustersPerEvent];
//  short tCluster_NumberOfHits[NMaxClustersPerEvent];
//  short tCluster_View[NMaxClustersPerEvent];
//  float tCluster_ChargeIntegral[NMaxHitsPerEvent];
//  short tCluster_StartChannel[NMaxHitsPerEvent];
//  short tCluster_StartTick[NMaxHitsPerEvent];
//  short tCluster_EndChannel[NMaxHitsPerEvent];
//  short tCluster_EndTick[NMaxHitsPerEvent];
//  float tCluster_StartCharge[NMaxHitsPerEvent];
//  float tCluster_StartAngle[NMaxHitsPerEvent];
//  float tCluster_EndCharge[NMaxHitsPerEvent];
//  float tCluster_EndAngle[NMaxHitsPerEvent];

  //Track variables
  short tNumberOfTracks;
//  short tTrackID[NMaxTracksPerEventTimesNViews];
  short tTrack_NumberOfHits[NMaxTracksPerEventTimesNViews];
  float tTrack_Length_Trajectory[NMaxTracksPerEventTimesNViews];
  float tTrack_Length_StraightLine[NMaxTracksPerEventTimesNViews];

  float tTrack_StartPoint_X[NMaxTracksPerEvent];
  float tTrack_StartPoint_Y[NMaxTracksPerEvent];
  float tTrack_StartPoint_Z[NMaxTracksPerEvent];
//  float tTrack_StartPoint_DistanceToBoundary[NMaxTracksPerEvent];
  float tTrack_EndPoint_X[NMaxTracksPerEvent];
  float tTrack_EndPoint_Y[NMaxTracksPerEvent];
  float tTrack_EndPoint_Z[NMaxTracksPerEvent];
//  float tTrack_EndPoint_DistanceToBoundary[NMaxTracksPerEvent];
  
  float tTrack_StartDirection_Theta[NMaxTracksPerEvent];
  float tTrack_StartDirection_Phi[NMaxTracksPerEvent];
//  float tTrack_StartDirection_X[NMaxTracksPerEvent];
//  float tTrack_StartDirection_Y[NMaxTracksPerEvent];
//  float tTrack_StartDirection_Z[NMaxTracksPerEvent];

  float tTrack_EndDirection_Theta[NMaxTracksPerEvent];
  float tTrack_EndDirection_Phi[NMaxTracksPerEvent];
//  float tTrack_EndDirection_X[NMaxTracksPerEvent];
//  float tTrack_EndDirection_Y[NMaxTracksPerEvent];
//  float tTrack_EndDirection_Z[NMaxTracksPerEvent];

//  float tTrack_PitchInViews[NMaxTracksPerEvent];
  short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2];

  //Track hit variables
  float tTrack_Hit_X[NMaxHitsPerEvent];
  float tTrack_Hit_Y[NMaxHitsPerEvent];
  float tTrack_Hit_Z[NMaxHitsPerEvent];
  float tTrack_Hit_Theta[NMaxHitsPerEvent];
  float tTrack_Hit_Phi[NMaxHitsPerEvent];
  float tTrack_dx_3DPosition[NMaxHitsPerEvent];
  float tTrack_ds_3DLocalTrackDirection[NMaxHitsPerEvent];
//  short tTrack_Hit_TPC[NMaxHitsPerEvent];
  short tTrack_Hit_View[NMaxHitsPerEvent];
//  short tTrack_Hit_Channel[NMaxHitsPerEvent];
//  float tTrack_Hit_PeakTime[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];
//  float tTrack_Hit_PeakHeight[NMaxHitsPerEvent];
//  float tTrack_Hit_StartTime[NMaxHitsPerEvent];
//  float tTrack_Hit_EndTime[NMaxHitsPerEvent];
//  float tTrack_Hit_Width[NMaxHitsPerEvent];
  float tTrack_Hit_GoodnessOfFit[NMaxHitsPerEvent];
//  short tTrack_Hit_Multiplicity[NMaxHitsPerEvent];

  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("Run",&tRun);
  rTree->SetBranchAddress("Subrun",&tSubrun);
  rTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  rTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
//  rTree->SetBranchAddress("EventTimeNanoseconds",&tEventTimeNanoseconds);
//  rTree->SetBranchAddress("IsData",&tIsData);

//  //Hit variables
  rTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);
//  rTree->SetBranchAddress("Hit_TPC",&tHit_TPC);
  rTree->SetBranchAddress("Hit_View",&tHit_View);
//  rTree->SetBranchAddress("Hit_Channel",&tHit_Channel);
//  rTree->SetBranchAddress("Hit_ChargeSummedADC",&tHit_ChargeSummedADC);
//  rTree->SetBranchAddress("Hit_ChargeIntegral",&tHit_ChargeIntegral);
//  rTree->SetBranchAddress("Hit_PeakHeight",&tHit_PeakHeight);
//  rTree->SetBranchAddress("Hit_PeakTime",&tHit_PeakTime);
//  rTree->SetBranchAddress("Hit_StartTime",&tHit_StartTime);
//  rTree->SetBranchAddress("Hit_EndTime",&tHit_EndTime);
//  rTree->SetBranchAddress("Hit_Width",&tHit_Width);
  rTree->SetBranchAddress("Hit_GoodnessOfFit",&tHit_GoodnessOfFit);
//  rTree->SetBranchAddress("Hit_Multiplicity",&tHit_Multiplicity);
  rTree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);
//  rTree->SetBranchAddress("Hit_ClusterID",&tHit_ClusterID);

//  //Cluster variables
//  rTree->SetBranchAddress("NumberOfClusters",&tNumberOfClusters);
//  rTree->SetBranchAddress("ClusterID",&tClusterID);
//  rTree->SetBranchAddress("Cluster_View",&tCluster_View);
//  rTree->SetBranchAddress("Cluster_NumberOfHits",&tCluster_NumberOfHits);
//  rTree->SetBranchAddress("Cluster_ChargeIntegral",&tCluster_ChargeIntegral);
//  rTree->SetBranchAddress("Cluster_StartChannel",&tCluster_StartChannel);
//  rTree->SetBranchAddress("Cluster_StartTick",&tCluster_StartTick);
//  rTree->SetBranchAddress("Cluster_EndChannel",&tCluster_EndChannel);
//  rTree->SetBranchAddress("Cluster_EndTick",&tCluster_EndTick);
//  rTree->SetBranchAddress("Cluster_StartCharge",&tCluster_StartCharge);
//  rTree->SetBranchAddress("Cluster_StartAngle",&tCluster_StartAngle);
//  rTree->SetBranchAddress("Cluster_EndCharge",&tCluster_EndCharge);
//  rTree->SetBranchAddress("Cluster_EndAngle",&tCluster_EndAngle);

//  //Track variables
  rTree->SetBranchAddress("NumberOfTracks",&tNumberOfTracks);
//  rTree->SetBranchAddress("TrackID",&tTrackID);
  rTree->SetBranchAddress("Track_NumberOfHits",&tTrack_NumberOfHits);
  rTree->SetBranchAddress("Track_Length_Trajectory",&tTrack_Length_Trajectory);
  rTree->SetBranchAddress("Track_Length_StraightLine",&tTrack_Length_StraightLine);

  rTree->SetBranchAddress("Track_StartPoint_X", &tTrack_StartPoint_X);
  rTree->SetBranchAddress("Track_StartPoint_Y", &tTrack_StartPoint_Y);
  rTree->SetBranchAddress("Track_StartPoint_Z", &tTrack_StartPoint_Z);
//  rTree->SetBranchAddress("Track_StartPoint_DistanceToBoundary", &tTrack_StartPoint_DistanceToBoundary);
  rTree->SetBranchAddress("Track_EndPoint_X", &tTrack_EndPoint_X);
  rTree->SetBranchAddress("Track_EndPoint_Y", &tTrack_EndPoint_Y);
  rTree->SetBranchAddress("Track_EndPoint_Z", &tTrack_EndPoint_Z);
//  rTree->SetBranchAddress("Track_EndPoint_DistanceToBoundary",&tTrack_EndPoint_DistanceToBoundary);

  rTree->SetBranchAddress("Track_StartDirection_Theta",&tTrack_StartDirection_Theta);
  rTree->SetBranchAddress("Track_StartDirection_Phi",&tTrack_StartDirection_Phi);
//  rTree->SetBranchAddress("Track_StartDirection_X", &tTrack_StartDirection_X);
//  rTree->SetBranchAddress("Track_StartDirection_Y", &tTrack_StartDirection_Y);
//  rTree->SetBranchAddress("Track_StartDirection_Z", &tTrack_StartDirection_Z);

  rTree->SetBranchAddress("Track_EndDirection_Theta",&tTrack_EndDirection_Theta);
  rTree->SetBranchAddress("Track_EndDirection_Phi",&tTrack_EndDirection_Phi);
//  rTree->SetBranchAddress("Track_EndDirection_X", &tTrack_EndDirection_X);
//  rTree->SetBranchAddress("Track_EndDirection_Y", &tTrack_EndDirection_Y);
//  rTree->SetBranchAddress("Track_EndDirection_Z", &tTrack_EndDirection_Z);

//  rTree->SetBranchAddress("Track_PitchInViews", &tTrack_PitchInViews);
  rTree->SetBranchAddress("Track_NumberOfHitsPerView",&tTrack_NumberOfHitsPerView);

//  //Track hit variables
  rTree->SetBranchAddress("Track_Hit_X", &tTrack_Hit_X);
  rTree->SetBranchAddress("Track_Hit_Y", &tTrack_Hit_Y);
  rTree->SetBranchAddress("Track_Hit_Z", &tTrack_Hit_Z);
  rTree->SetBranchAddress("Track_Hit_LocalTrackDirection_Theta", &tTrack_Hit_Theta);
  rTree->SetBranchAddress("Track_Hit_LocalTrackDirection_Phi", &tTrack_Hit_Phi);
  rTree->SetBranchAddress("Track_Hit_ds_3DPosition", &tTrack_dx_3DPosition);
  rTree->SetBranchAddress("Track_Hit_ds_LocalTrackDirection", &tTrack_ds_3DLocalTrackDirection);
//  rTree->SetBranchAddress("Track_Hit_TPC", &tTrack_Hit_TPC);
  rTree->SetBranchAddress("Track_Hit_View", &tTrack_Hit_View);
//  rTree->SetBranchAddress("Track_Hit_Channel", &tTrack_Hit_Channel);
//  rTree->SetBranchAddress("Track_Hit_PeakTime", &tTrack_Hit_PeakTime);
  rTree->SetBranchAddress("Track_Hit_ChargeSummedADC", &tTrack_Hit_ChargeSummedADC);
  rTree->SetBranchAddress("Track_Hit_ChargeIntegral", &tTrack_Hit_ChargeIntegral);
//  rTree->SetBranchAddress("Track_Hit_PeakHeight", &tTrack_Hit_PeakHeight);
//  rTree->SetBranchAddress("Track_Hit_StartTime", &tTrack_Hit_StartTime);
//  rTree->SetBranchAddress("Track_Hit_EndTime", &tTrack_Hit_EndTime);
//  rTree->SetBranchAddress("Track_Hit_Width", &tTrack_Hit_Width);
  rTree->SetBranchAddress("Track_Hit_GoodnessOfFit", &tTrack_Hit_GoodnessOfFit);
//  rTree->SetBranchAddress("Track_Hit_Multiplicity", &tTrack_Hit_Multiplicity);
  if( to_read == 0){
    to_read = NEntries;
  }
  double Efield = -1.;
  
  for(int i=0; i<to_read; i++) //Event loop
  {
    rTree->GetEntry(i);
    if(i==0){
      tstart = tEventTimeSeconds;
      TMyFileHeader header = load_run_header(tRun);
      if(header.GetRun() == -1){return;}
      Efield = header.GetAmplification();
    }
    if(i==to_read-1){
      tend = tEventTimeSeconds;
    }

    //initialize classes (not pointers for the moment)

    int a=0; //Need this counter to remember where we are in this array.
    for(int j=0; j<tNumberOfTracks; j++){ //Track loop
    
      track dummy_track;
      dummy_track.nhitsview[0] = 0; dummy_track.nhitsview[1] = 0;
      dummy_track.throughgoing = false;
      dummy_track.throughgoingX = false;
      dummy_track.throughgoingY = false;
      dummy_track.throughgoingZ = false;
      dummy_track.on_edge = false;
      dummy_track.on_edgeX = false;
      dummy_track.on_edgeY = false;
      dummy_track.on_edgeZ = false;
      dummy_track.highway_selected = false;
      dummy_track.blc = false;

      dummy_track.run = tRun;
      dummy_track.subrun = tSubrun;
      dummy_track.event = tEventNumberInRun;
      dummy_track.event_ntracks = tNumberOfTracks;
      dummy_track.id = j;
      dummy_track.start_x = tTrack_StartPoint_X[j];
      dummy_track.start_y = tTrack_StartPoint_Y[j];
      dummy_track.start_z = tTrack_StartPoint_Z[j];

      dummy_track.end_x = tTrack_EndPoint_X[j];
      dummy_track.end_y = tTrack_EndPoint_Y[j];
      dummy_track.end_z = tTrack_EndPoint_Z[j];
      
      
      if( max(dummy_track.end_x, dummy_track.start_x) > MaxX and min(dummy_track.end_x, dummy_track.start_x) < MinX ){
        dummy_track.throughgoingX = true;
      } //end if x
      if( max(dummy_track.end_x, dummy_track.start_x) > MaxX or min(dummy_track.end_x, dummy_track.start_x) < MinX ){
        dummy_track.on_edgeX = true;
      }
      if( max(dummy_track.end_y, dummy_track.start_y) > MaxY and min(dummy_track.end_y, dummy_track.start_y) < MinY ){
        dummy_track.throughgoingY = true;
      } //end if z
      if( max(dummy_track.end_y, dummy_track.start_y) > MaxY or min(dummy_track.end_y, dummy_track.start_y) < MinY ){
        dummy_track.on_edgeY = true;
      }
      if( max(dummy_track.end_z, dummy_track.start_z) > MaxZ and min(dummy_track.end_z, dummy_track.start_z) < MinZ ){
        dummy_track.throughgoingZ = true;
      } //end if y
      if( max(dummy_track.end_z, dummy_track.start_z) > MaxZ or min(dummy_track.end_z, dummy_track.start_z) < MinZ ){
        dummy_track.on_edgeZ = true;
      }
      if(dummy_track.throughgoingX or dummy_track.throughgoingY or dummy_track.throughgoingZ){
        dummy_track.throughgoing = true;
      }
      if(dummy_track.on_edgeX or dummy_track.on_edgeY or dummy_track.on_edgeZ){
        dummy_track.on_edge = true;
      }
      
      int lem0 = find_lem(0, dummy_track.start_z);
      int lem1 = find_lem(0, dummy_track.end_z);
      if(isGood_lem(lem0) and isGood_lem(lem1) ){dummy_track.blc = true;}
      if(highway){
        if(tracks_selected_by_highway.size() == 0){
          cout << "ERROR in Read_tree_June : please load highway output" << endl;
          return;
        }
        int track_ID = 1000*(1+dummy_track.event) + dummy_track.id;
        if( find(tracks_selected_by_highway.begin(), tracks_selected_by_highway.end(), track_ID) != tracks_selected_by_highway.end() ){
          dummy_track.highway_selected = true;
        }
      }
      
      dummy_track.start_theta = fabs(tTrack_StartDirection_Theta[j]);
      if( dummy_track.start_theta > 180 ){
        dummy_track.start_theta = 360 - dummy_track.start_theta;
      }
      if( dummy_track.start_theta < 90 ){
        dummy_track.start_theta = 180 - dummy_track.start_theta;
      }
      dummy_track.start_phi = tTrack_StartDirection_Phi[j];
//      dummy_track.start_vx = tTrack_StartDirection_X[j];
//      dummy_track.start_vy = tTrack_StartDirection_Y[j];
//      dummy_track.start_vz = tTrack_StartDirection_Z[j];
      
      dummy_track.end_theta = fabs(tTrack_EndDirection_Theta[j]);
      if( dummy_track.end_theta > 180 ){
        dummy_track.end_theta = 360 - dummy_track.end_theta;
      }
      if( dummy_track.end_theta < 90 ){
        dummy_track.end_theta = 180 - dummy_track.end_theta;
      }
      dummy_track.end_phi = tTrack_EndDirection_Phi[j];
//      dummy_track.end_vx = tTrack_EndDirection_X[j];
//      dummy_track.end_vy = tTrack_EndDirection_Y[j];
//      dummy_track.end_vz = tTrack_EndDirection_Z[j];

      dummy_track.length = tTrack_Length_StraightLine[j];

      dummy_track.theta = get_theta(dummy_track);
      dummy_track.phi = get_phi(dummy_track);
      dummy_track.long_in_view = -1;
      if(get_reduced_phi(dummy_track) < 20 ){
        dummy_track.long_in_view = 0;
      }
      if(get_reduced_phi(dummy_track) > 70 ){
        dummy_track.long_in_view = 1;
      }
      
      for(int l=a; l < (a+tTrack_NumberOfHits[j]); l++){
        hit dummy_hits;
        double myds = 0;
        double mytheta = fabs(tTrack_Hit_Theta[l]);
        double myphi = tTrack_Hit_Phi[l];
        if(mytheta > 180){mytheta = 360-mytheta;}
        if(mytheta < 90){mytheta = 180-mytheta;}
        if(tTrack_Hit_View[l] == 0){
          myds = pitch / (TMath::Sin(mytheta*TMath::Pi()/180.)*TMath::Sin(myphi*TMath::Pi()/180.));
        }
        if(tTrack_Hit_View[l] == 1){
          myds = pitch / (TMath::Sin(mytheta*TMath::Pi()/180.)*TMath::Cos(myphi*TMath::Pi()/180.));
        }
        
        dummy_hits.run = tRun;
        dummy_hits.subrun = tSubrun;
        dummy_hits.event = tEventNumberInRun;
        dummy_hits.view = tTrack_Hit_View[l];
        dummy_hits.track_id = dummy_track.id;
        dummy_hits.dq_sum = tTrack_Hit_ChargeSummedADC[l]/ADC2CHARGE[tTrack_Hit_View[l]];
        dummy_hits.dq_integral = tTrack_Hit_ChargeIntegral[l]/ADC2CHARGE[tTrack_Hit_View[l]];
        dummy_hits.dx_3D = tTrack_dx_3DPosition[l];
        dummy_hits.dx_local = tTrack_ds_3DLocalTrackDirection[l];
        dummy_hits.sp_x = tTrack_Hit_X[l];
        dummy_hits.sp_y = tTrack_Hit_Y[l];
        dummy_hits.sp_z = tTrack_Hit_Z[l];
        dummy_hits.sp_theta = mytheta;
        dummy_hits.sp_phi = myphi;
        dummy_hits.dx_phil = myds;
        dummy_hits.purity_correction = correct_dx_for_lifetime(50-tTrack_Hit_X[l]);
        dummy_hits.lem  = find_lem(dummy_hits.sp_y, dummy_hits.sp_z);
        dummy_hits.gain_density_correction_factor = gain_correction_for_rho(get_density_for_hit(tEventTimeSeconds, tRun), Efield);
        dummy_hits.GoF = tHit_GoodnessOfFit[l];
        dummy_track.nhitsview[dummy_hits.view]++;

        dummy_track.hits_trk.push_back(dummy_hits);
      } //end hits
      tracks.push_back(dummy_track);
      a+=tTrack_NumberOfHits[j];
    }//end tracks
  }//Event loop
  
  
  return;
}

//Analysis *********************************************************************

void FWHM(double &st, double &end, TH1D *h){
  //find fit paramter using FWHM. (same approach as uBoone)
  st = 0;
  end = 50;
  int bin_min = 3;
  double x_max = find_max_bin( h, bin_min ); // x position of max of histo
  double bin_max = h->FindBin( x_max ); // bin number of xmax
  double n_max = h->GetBinContent( bin_max ); // max of histo
  double th=0.5*n_max;
  for(int bin=bin_min; bin<bin_max; bin++ ){//loop from first bin to bin of max of histo
    if( th < h->GetBinContent(bin) or bin_max-bin < 3){//set start x to correspond to half of histo maximum, at least 3 bins away from maximum
      st=h->GetBinCenter( bin-1 )+h->GetBinWidth(bin-1)/2.;
      break;
    }
  }
  for(int bin=bin_max; bin<h->GetNbinsX(); bin++ ){
    if( th > h->GetBinContent(bin) and bin-bin_max > 3){//set end x to correspond to half of histo maximum, at least 4 bins away from maximum
      end=h->GetBinCenter( bin+1 )+h->GetBinWidth(bin-1)/2.;
      break;
    }
  }
  //adjustements if start and end are nonsense
  if ( st <= 0 and x_max > 4 ){
    st = x_max-4;
  }
  if ( x_max - st < 4 and x_max > 4){
    st = x_max-4;
  }
  if ( end <= 0 ){
    end = x_max+5;
  }
  if ( end - x_max < 5 ){
    end = x_max+5;
  }
  return;
}

double find_closest_fit_value_in_bin(TH1D* myhisto, int bin, TF1 *myfunc){
  double low_x = myhisto->GetBinCenter(bin)-myhisto->GetBinWidth(bin)/2.;
  double up_x = myhisto->GetBinCenter(bin)+myhisto->GetBinWidth(bin)/2.;
  double biny = myhisto->GetBinContent(bin);
  double low_y = myfunc->Eval(low_x);
  double up_y = myfunc->Eval(up_x);
  double best_value = -1;
       if(low_y > biny and low_y < up_y){best_value = low_y;}
  else if(up_y  < biny and low_y < up_y){best_value = up_y;}
  else if(low_y < biny and low_y > up_y){best_value = low_y;}
  else if(up_y  > biny and low_y > up_y){best_value = up_y;}
  else if(low_y == low_x){best_value = low_y;}
  else {best_value = biny;}
  return best_value;
}

double homemade_get_pvalue(TF1 *myfunc, TH1D *myhisto){
  double startx, endx, dx;
  myfunc->GetRange(startx,endx);
  int startbin = myhisto->FindBin(startx);
  int endbin = myhisto->FindBin(endx);
  int binrange = endbin-startbin;
  double chisquare = 0;
  for(int i = startbin; i < endbin; i++){
    double yhisto = myhisto->GetBinContent(i);
    double yfunc = find_closest_fit_value_in_bin(myhisto, i, myfunc);
    chisquare += pow(yhisto-yfunc,2)/(yhisto);
  }
  return TMath::Prob(chisquare, binrange-1);
}

double truncated_mean(std::vector<double> &vec, double rmlow, double rmhigh ){
  //**************************************************************************//
  // truncated mean
  //   rmlow  - fraction of events to remove at the low part of the distribution
  //   rmhigh - fraction of events to remove at the high part of the distribution
  //

 if(vec.empty()) return -999;

 std::vector<double> vtmp(vec);

 // sort it
 std::sort(vtmp.begin(), vtmp.end());

 size_t n = vtmp.size();
 size_t nlow, nhigh;

  //evaluate the number of entries to trim from the start and from the end
 // TODO: proceeding in this way, rounding the entry, we might be biased
 nlow  = (size_t)(rmlow * n);
 nhigh = n - (size_t)(rmhigh * n);

 // some basic checks
 if(rmlow <= 0 or rmlow >= 1) nlow = 0;
 if(rmhigh <= 0 or rmhigh >= 1) nhigh = n;

 if( (rmlow + rmhigh) >= 0.99)
   {
     cout<<"ERROR: truncated_mean() the fractions to remove are to high"<<endl;
     nlow  = 0;
     nhigh = n;
   }

 //cout<<“Subsample “<<nlow<<“, “<<nhigh<<” out of “<<n<<endl;

 // should not happen though
 if( nlow >= n ) nlow = 0;
 if( nhigh >= n ) nhigh = n;


 double trmean = 0.0;

 for(size_t i=nlow;i<nhigh;i++)
   trmean += vtmp[i];

 trmean /= ((float)nhigh - (float)nlow);
 return trmean;
}

double find_max_bin(TH1D *hist, int min_bin){
  //return the position of the bin with the largest number of entries
  int bin_max =0; double n_bin_max =0; double x_max =0.;

  for(int bin = min_bin; bin < (hist->GetSize()-2); bin++ ){
    if(hist->GetBinContent(bin) > n_bin_max){
      n_bin_max = hist->GetBinContent(bin);
      bin_max = bin;
    } //end if
  } //end for bin

  x_max = hist->GetXaxis()->GetBinCenter(bin_max);
  return x_max;
}

double langaufun(double *x, double *par){

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;

      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1D *his, double *fitrange, vector<double> startvalues,
 vector<double> parlimitslo, vector<double> parlimitshi, vector<double> &fitparams, vector<double> &fiterrors, double &pvalue, double pvaluelim, bool find_best, bool gauss){
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf
   //   If find_best == true, is returned the fit with the lower chisquare within multiple ranges.
  
  double mpshift  = -0.22278298;
  TString FunName = TString::Format("Fitfcn_%s",his->GetName());
  vector<string> parname;
  parname.push_back("Width"); parname.push_back("MP"); parname.push_back("Area"); parname.push_back("GSigma");

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName.Data());
  if (ffitold) delete ffitold;
  TF1 *ffit;
  if(!gauss){
    ffit = new TF1(FunName.Data(), langaufun, fitrange[0], fitrange[1],4);
  }
  else{
    double mpc = startvalues[1] - mpshift * startvalues[0];
    startvalues[1] = mpc;
    ffit = new TF1(FunName.Data(), "gaus", fitrange[0], fitrange[1]);
    startvalues.pop_back();
    parlimitslo.pop_back();
    parlimitshi.pop_back();
    parname.pop_back();
  }

  for( unsigned int i=0; i<startvalues.size(); i++ ){
    ffit->SetParameter(i, startvalues[i]);
    ffit->SetParName(i, parname[i].data());
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  //lower bound varies from 0.0*mean to 0.5*mean in septs of 0.05
  //upper bound varies from 0.9*mean to 3.1*mean in
  double max_p = 0; /*double min_chi =9999.;*/ double min_j =0.; double min_i =0.;
  double x_max = find_max_bin(his);

  if( find_best ){
    for( double i=0.1; i < fitrange[1]; i+=0.05 ){
      for( double j=fitrange[0]; j<3.0; j+=0.15 ){
        ffit->SetRange( i*his->GetMean(),j*his->GetMean() );
        his->Fit(FunName.Data(),"RBSQ");   // fit within specified range, use ParLimits, do not plot
        if( i*his->GetMean() < x_max and j*his->GetMean() > x_max ){
//          if( ( ffit->GetChisquare()/ffit->GetNDF() ) < min_chi ){
//            min_chi = ffit->GetChisquare()/ffit->GetNDF();
//            min_j = j; min_i =i;
//          }
          if( ffit->GetProb() > max_p ){
            max_p = ffit->GetProb();
            min_j = j; min_i =i;
          }
        }
      }//end j
    }//end i
    ffit->SetRange( min_i*his->GetMean(),min_j*his->GetMean() );
    his->Fit(FunName.Data(),"LRBSQ");
  }//end find_best
  else{
    his->Fit(FunName.Data(),"LRBSQ");
  }
  pvalue = ffit->GetProb();
//  pvalue = homemade_get_pvalue(ffit, his);
//  if(pvalue < pvaluelim){
//    double bin_max = his->FindBin( find_max_bin( his, 4 ) );
//    fitrange[0] = his->GetBinCenter(bin_max-4);
//    fitrange[1] = his->GetBinCenter(bin_max+7);
//    ffit->SetRange(fitrange[0], fitrange[1]);
//    his->GetListOfFunctions()->Clear();
//    his->Fit(FunName.Data(),"LRBSQ");
//    pvalue = homemade_get_pvalue(ffit, his);
//  }
  for (unsigned int i=0; i<startvalues.size(); i++) {
    fitparams.push_back(ffit->GetParameter(i));    // obtain fit parameters
    fiterrors.push_back(ffit->GetParError(i));     // obtain fit parameter errors
  }

  return (ffit);              // return fit function
}

void plot_tracks(vector<track> tracks, string cut_type, int maxtracks, double dqds_cut){
  
  if(tracks.size() == 0){
    return;
  }
  
  int run = tracks[0].run;
  int trk = 0;
  
  string outfile = SelectTrack_Output + cut_type + "/plots/" + to_string(run);
  check_and_mkdir(outfile);
  outfile = outfile + "/tracks";
  if(maxtracks > 0){
    outfile = outfile + "_0_to_" + to_string(maxtracks);
  }
  if(dqds_cut != 0){
    if(dqds_cut < 0){
      outfile = outfile + "_dQds_inf_to_" + to_string_with_precision(dqds_cut,2);
    }
    else{
      outfile = outfile + "_dQds_sup_to_" + to_string_with_precision(dqds_cut,2);
    }
  }
  TH1D dQds0("dQds0","view 0;dQds (fC/cm)",100,0,50);
  dQds0.GetXaxis()->CenterTitle();
  dQds0.GetYaxis()->CenterTitle();
  TH1D dQds1("dQds1","view 1;dQds (fC/cm)",100,0,50);
  dQds1.GetXaxis()->CenterTitle();
  dQds1.GetYaxis()->CenterTitle();
  TH2D XY("XY_front_view","Front (XY) view;y (cm);x (cm);dQds (fC/cm)",100/pitch,-50,50,100/pitch,-50,50);
  XY.GetXaxis()->CenterTitle();
  XY.GetYaxis()->CenterTitle();
  TH2D XZ("XZ_side_view","Side (XZ) view;z (cm);x (cm);dQds (fC/cm)",300/pitch,0,300,100/pitch,-50,50);
  XZ.GetXaxis()->SetLabelSize(0.07); 
  XZ.GetXaxis()->SetTitleSize(0.07);
  XZ.GetXaxis()->CenterTitle();
  XZ.GetYaxis()->CenterTitle();
  XZ.GetYaxis()->SetLabelSize(0.07); 
  XZ.GetYaxis()->SetTitleSize(0.07);
  TH2D YZ("YZ_top_view","Top (YZ) view;z (cm);y (cm);dQds (fC/cm)",300/pitch,0,300,100/pitch,-50,50);
  YZ.GetXaxis()->SetLabelSize(0.07); 
  YZ.GetXaxis()->SetTitleSize(0.07);
  YZ.GetXaxis()->CenterTitle();
  YZ.GetYaxis()->CenterTitle();
  YZ.GetYaxis()->SetLabelSize(0.07); 
  YZ.GetYaxis()->SetTitleSize(0.07);
  TH3D XYZ("XYZ","XYZ;z (cm);y (cm);x (cm)",300/pitch,0,300,100/pitch,-50,50,100/pitch,-50,50);
  XYZ.GetXaxis()->CenterTitle();
  XYZ.GetYaxis()->CenterTitle();
  XYZ.GetZaxis()->CenterTitle();
  XY.SetStats(false);
  XZ.SetStats(false);
  YZ.SetStats(false);
  XYZ.SetStats(false);
  
  for( auto atrack : tracks){
  
    trk++;
    if(trk == maxtracks){
      break;
    }
    #if verbose
    cout << "Plotting track " << trk << endl;
    #endif 
    TGraph onetheta0, onetheta1, onephi0, onephi1;
    TH1D onedQ0(string("dQ0_Event_" + to_string(atrack.event) + "_track_" + to_string(atrack.id)).data(),string("view 0. Total charge: " + to_string((int)atrack.sumQ[0]) + "fC;dQ (fC)").data(),100,0,200);
    onedQ0.GetXaxis()->CenterTitle();
    onedQ0.GetYaxis()->CenterTitle();
    TH1D onedQ1(string("dQ1_Event_" + to_string(atrack.event) + "_track_" + to_string(atrack.id)).data(),string("view 1. Total charge: " + to_string((int)atrack.sumQ[1]) + "fC;dQ (fC)").data(),100,0,200);
    onedQ1.GetXaxis()->CenterTitle();
    onedQ1.GetYaxis()->CenterTitle();
    TH1D oneds0(string("ds0_Event_" + to_string(atrack.event) + "_track_" + to_string(atrack.id)).data(),"view 0;ds (cm)",100,0,10);
    oneds0.GetXaxis()->CenterTitle();
    oneds0.GetYaxis()->CenterTitle();
    TH1D oneds1(string("ds1_Event_" + to_string(atrack.event) + "_track_" + to_string(atrack.id)).data(),"view 1;ds (cm)",100,0,10);
    oneds1.GetXaxis()->CenterTitle();
    oneds1.GetYaxis()->CenterTitle();
    TH1D onedQds0(string("dQds0_Event_" + to_string(atrack.event) + "_track_" + to_string(atrack.id)).data(),"view 0;dQds (fC/cm)",100,0,50);
    onedQds0.GetXaxis()->CenterTitle();
    onedQds0.GetYaxis()->CenterTitle();
    TH1D onedQds1(string("dQds1_Event_" + to_string(atrack.event) + "_track_" + to_string(atrack.id)).data(),"view 1;dQds (fC/cm)",100,0,50);
    onedQds1.GetXaxis()->CenterTitle();
    onedQds1.GetYaxis()->CenterTitle();
    TH2D oneXY("oneXY_front_view","Front (XY) view;y (cm);x (cm);dQds (fC/cm)",500,-50,50,500,-50,50);
    oneXY.GetXaxis()->CenterTitle();
    oneXY.GetYaxis()->CenterTitle();
    TH2D oneXZ("oneXZ_side_view","Side (XZ) view;z (cm);x (cm);dQds (fC/cm)",1500,0,300,500,-50,50);
    oneXZ.GetXaxis()->CenterTitle();
    oneXZ.GetXaxis()->SetLabelSize(0.07); 
    oneXZ.GetXaxis()->SetTitleSize(0.07);
    oneXZ.GetYaxis()->CenterTitle();
    oneXZ.GetYaxis()->SetLabelSize(0.07); 
    oneXZ.GetYaxis()->SetTitleSize(0.07);
    TH2D oneYZ("oneYZ_top_view","Top (YZ) view;z (cm);y (cm);dQds (fC/cm)",1500,0,300,500,-50,50);
    oneYZ.GetXaxis()->SetLabelSize(0.07); 
    oneYZ.GetXaxis()->SetTitleSize(0.07);
    oneYZ.GetXaxis()->CenterTitle();
    oneYZ.GetYaxis()->CenterTitle();
    oneYZ.GetYaxis()->SetLabelSize(0.07); 
    oneYZ.GetYaxis()->SetTitleSize(0.07);
    TH3D oneXYZ("oneXYZ","XYZ;z (cm);y (cm);x (cm)",600,0,300,200,-50,50,200,-50,50);
    oneXYZ.GetXaxis()->CenterTitle();
    oneXYZ.GetYaxis()->CenterTitle();
    oneXYZ.SetStats(false);
    oneXYZ.SetStats(false);
    oneXYZ.SetStats(false);
  
    for(auto ahit : atrack.hits_trk){
      double dq;
      double ds;
      if(method_ds == "3D"){ds = ahit.dx_3D;}
      else if(method_ds == "local"){ds = ahit.dx_local;}
      else{ds = ahit.dx_phil;}
      if(method_dQ == "sum"){dq = ahit.dq_sum;}
      else{dq = ahit.dq_integral;}
      double dqds = dq/ds;
      bool Continue = false;
      if(dqds_cut != 0){
        if(dqds_cut < 0){
          if(dqds > fabs(dqds_cut)){Continue = true;}
        }
        if(dqds_cut > 0){
          if(dqds < dqds_cut){Continue = true;}
        }
      }
      if(Continue){continue;}
      if(ahit.view == 0){
        onedQ0.Fill(dq);
        oneds0.Fill(ds);
        onedQds0.Fill(dqds);
        dQds0.Fill(dqds);
        onetheta0.SetPoint(onetheta0.GetN(),onetheta0.GetN(),ahit.sp_theta);
        onephi0.SetPoint(onephi0.GetN(),onephi0.GetN(),ahit.sp_phi);
      }
      else{
        onedQ1.Fill(dq);
        oneds1.Fill(ds);
        onedQds1.Fill(dqds);
        dQds1.Fill(dqds);
        onetheta1.SetPoint(onetheta1.GetN(),onetheta1.GetN(),ahit.sp_theta);
        onephi1.SetPoint(onephi1.GetN(),onephi1.GetN(),fabs(ahit.sp_phi));
      }
      oneXY.Fill(ahit.sp_y,ahit.sp_x, dqds);
      oneXZ.Fill(ahit.sp_z,ahit.sp_x, dqds);
      oneYZ.Fill(ahit.sp_z,ahit.sp_y, dqds);
      oneXYZ.Fill(ahit.sp_z,ahit.sp_y,ahit.sp_x, dqds);
      XY.Fill(ahit.sp_y,ahit.sp_x, dqds);
      XZ.Fill(ahit.sp_z,ahit.sp_x, dqds);
      YZ.Fill(ahit.sp_z,ahit.sp_y, dqds);
      XYZ.Fill(ahit.sp_z,ahit.sp_y,ahit.sp_x, dqds);
    }
    TCanvas onetrack_can("tracks","tracks",2000,3000);
      onetrack_can.Divide(1,6);
      onetrack_can.cd(1);
      TPad *oneXY_pad = new TPad("oneXY_pad","oneXY_pad",0,0,0.25,1);
      oneXY_pad->Draw();
      oneXY_pad->cd();
        oneXY.Draw("COLZ");
      onetrack_can.cd(1);
      TPad *oneXZ_pad = new TPad("oneXZ_pad","oneXZ_pad",0.25,0,1,0.5);
      oneXZ_pad->Draw();
      oneXZ_pad->cd();
        oneXZ.Draw("COLZ");
      onetrack_can.cd(1);
      TPad *oneYZ_pad = new TPad("oneYZ_pad","oneYZ_pad",0.25,0.51,1,1);
      oneYZ_pad->Draw();
      oneYZ_pad->cd();
        oneYZ.Draw("COLZ");
      onetrack_can.cd(2);
      TPad *onedQ0_pad = new TPad("onedQ0_pad","onedQ0_pad",0,0,0.5,1);
      onedQ0_pad->Draw();
      onedQ0_pad->cd();
        onedQ0.Draw();
      onetrack_can.cd(2);
      TPad *onedQ1_pad = new TPad("onedQ1_pad","onedQ1_pad",0.5,0,1,1);
      onedQ1_pad->Draw();
      onedQ1_pad->cd();
        onedQ1.Draw();
      onetrack_can.cd(3);
      TPad *oneds0_pad = new TPad("oneds0_pad","oneds0_pad",0,0,0.5,1);
      oneds0_pad->Draw();
      oneds0_pad->cd();
        oneds0.Draw();
      onetrack_can.cd(3);
      TPad *oneds1_pad = new TPad("oneds1_pad","oneds1_pad",0.5,0,1,1);
      oneds1_pad->Draw();
      oneds1_pad->cd();
        oneds1.Draw();
      onetrack_can.cd(4);
      TPad *onedQds0_pad = new TPad("onedQds0_pad","onedQds0_pad",0,0,0.5,1);
      onedQds0_pad->Draw();
      onedQds0_pad->cd();
        onedQds0.Draw();
      onetrack_can.cd(4);
      TPad *onedQds1_pad = new TPad("onedQds1_pad","onedQds1_pad",0.5,0,1,1);
      onedQds1_pad->Draw();
      onedQds1_pad->cd();
        onedQds1.Draw();
      onetrack_can.cd(5);
      TPad *onetheta0_pad = new TPad("onetheta0_pad","onetheta0_pad",0,0,0.5,1);
      onetheta0_pad->Draw();
      onetheta0_pad->cd();
        onetheta0.SetTitle("view 0;hit number;theta (rad)");
        onetheta0.Draw("AP");
        onetheta0.SetMinimum(90);onetheta0.SetMaximum(180);
        onetheta0.Draw("AP");
        onetheta0.GetXaxis()->CenterTitle();
        onetheta0.GetYaxis()->CenterTitle();
      onetrack_can.cd(5);
      TPad *onetheta1_pad = new TPad("onetheta1_pad","onetheta1_pad",0.5,0,1,1);
      onetheta1_pad->Draw();
      onetheta1_pad->cd();
        onetheta1.SetTitle("view 1;hit number;theta (rad)");
        onetheta1.Draw("AP");
        onetheta1.SetMinimum(90);onetheta1.SetMaximum(180);
        onetheta1.Draw("AP");
        onetheta1.GetXaxis()->CenterTitle();
        onetheta1.GetYaxis()->CenterTitle();
      onetrack_can.cd(6);
      TPad *onephi0_pad = new TPad("onephi0_pad","onephi0_pad",0,0,0.5,1);
      onephi0_pad->Draw();
      onephi0_pad->cd();
        onephi0.SetTitle("view 0;hit number;phi (rad)");
        onephi0.Draw("AP");
        onephi0.SetMinimum(-180);onephi0.SetMaximum(180);
        onephi0.Draw("AP");
        onephi0.GetXaxis()->CenterTitle();
        onephi0.GetYaxis()->CenterTitle();
      onetrack_can.cd(6);
      TPad *onephi1_pad = new TPad("onephi1_pad","onephi1_pad",0.5,0,1,1);
      onephi1_pad->Draw();
      onephi1_pad->cd();
        onephi1.SetTitle("view 1;hit number;phi (rad)");
        onephi1.Draw("AP");
        onephi1.SetMinimum(-180);onephi1.SetMaximum(180);
        onephi1.Draw("AP");
        onephi1.GetXaxis()->CenterTitle();
        onephi1.GetYaxis()->CenterTitle();
    if(trk == 1){
      onetrack_can.Print(string(outfile+".pdf(").data(),"pdf");
    }
    else{
      onetrack_can.Print(string(outfile+".pdf").data(),"pdf");
    }
  }
  
//  string outfileroot = outfile + ".root";
//  TFile ofile(outfileroot.data(), "RECREATE");
  TCanvas track_can("tracks","tracks",2000,1000);
    TPad *XY_pad = new TPad("XY_pad","XY_pad",0,0.5,0.25,1);
    XY_pad->Draw();
    XY_pad->cd();
      XY.Draw("COLZ");
    track_can.cd();
    TPad *XZ_pad = new TPad("XZ_pad","XZ_pad",0.25,0.51,1,0.75);
    XZ_pad->Draw();
    XZ_pad->cd();
      XZ.Draw("COLZ");
    track_can.cd();
    TPad *YZ_pad = new TPad("YZ_pad","YZ_pad",0.25,0.76,1,1);
    YZ_pad->Draw();
    YZ_pad->cd();
      YZ.Draw("COLZ");
    track_can.cd();
    TPad *dQds0_pad = new TPad("dQds0_pad","dQds0_pad",0,0,0.5,0.5);
    dQds0_pad->Draw();
    dQds0_pad->cd();
      dQds0.Draw();
    track_can.cd();
    TPad *dQds1_pad = new TPad("dQds1_pad","dQds1_pad",0.5,0,1,0.5);
    dQds1_pad->Draw();
    dQds1_pad->cd();
      dQds1.Draw();
//  track_can.Write();
  track_can.Print(string(outfile+".pdf").data(),"pdf");
//  XY.Write();
//  XZ.Write();
//  XZ.Write();
//  dQds0.Write();
//  dQds1.Write();
  
  
  TCanvas track_can_XYZ("tracks_3Dview","tracks 3D view",1500,500);
  track_can_XYZ.cd();
    XYZ.Draw("COLZ");
//  track_can_XYZ.Write();
  track_can_XYZ.Print(string(outfile+".pdf)").data(),"pdf");
//  XYZ.Write();
  
//  ofile.Close();
  return;
}

//////////not 311-related/////////////

inline vector<string> glob(const string& pat){
  glob_t glob_result;
  glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> ret;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
      string full_path = string(glob_result.gl_pathv[i]);
      string basename = full_path.substr(full_path.find_last_of("/")+1);
      ret.push_back(basename.data());
  }
  globfree(&glob_result);
  return ret;
}

inline bool ExistTest (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

vector<double> fit_dQds(TH1D *hdQds, bool gauss, int min_number_of_hits, double pvaluelim, double sigmalim, TGraphErrors* graph, float x, float xer){
  TF1 *function;
  if(hdQds->GetEntries() < min_number_of_hits){
    #if verbose
    cout << "    Ignoring " << hdQds->GetName() << ": too few hits (" << hdQds->GetEntries() << ")" << endl;
    #endif
    return {-1,-1};
  }

  //do fit
  #if verbose
  cout << "    Fit: " << hdQds->GetName() << endl;
  #endif

  double fr[2];
  vector<double> sv, pllo, plhi, fp, fpe;
  double chisqr, pvalue;

  //initialize fit quantities
  FWHM(fr[0], fr[1], hdQds);
  pllo.push_back(0.1*hdQds->GetStdDev());
  pllo.push_back(0.1*hdQds->GetMean());
  pllo.push_back(0.1*hdQds->Integral());
  pllo.push_back(0.1*hdQds->GetStdDev());
  plhi.push_back(10*hdQds->GetStdDev());
  plhi.push_back(10*hdQds->GetMean());
  plhi.push_back(100*hdQds->Integral());
  plhi.push_back(10*hdQds->GetStdDev());
  sv.push_back(hdQds->GetStdDev());
  sv.push_back(hdQds->GetMean());
  sv.push_back(hdQds->Integral());
  sv.push_back(hdQds->GetStdDev());
  function = langaufit(hdQds,fr,sv,pllo,plhi,fp,fpe,pvalue,pvaluelim, false, gauss);
  if( pvalue < pvaluelim or pvalue != pvalue ){
    #if verbose
    cout << "    Bad fit for " << hdQds->GetName() << "(pvalue=" << pvalue << ")" << endl;
    #endif
    return {-1,-1};
  }
  else if(fpe[1]/mpv_cosmics > sigmalim){
    #if verbose
    cout << "    Bad fit for " << hdQds->GetName() << "(sigma=" << fpe[1]/mpv_cosmics << ")" << endl;
    #endif
    return {-1,-1};
  }
//  double MPV = function->GetMaximumX();
//  double xmax, xmin;
//  function->GetRange(xmin,xmax);
//  int nstep = 0;
//  while ( (MPV == xmin or MPV == xmax) and nstep < 5 ){
//    function->SetRange(xmin-5,xmax+5);
//    function->GetRange(xmin,xmax);
//    MPV = function->GetMaximumX();
//    nstep++;
//  }
//  if(nstep == 5){
//    #if verbose
//    cout << "    Bad fit for " << hdQds->GetName() << " (MPV is equal to one of fit function ranges)" << endl;
//    #endif
//    return {-1,-1};
//  }
  ((TF1*)hdQds->GetListOfFunctions()->At(0))->SetLineColor(600);
  if( graph != 0 ){
    //insert points in graph
    graph->SetPoint(graph->GetN(), x, fp[1]/(corr*mpv_cosmics));
    graph->SetPointError(graph->GetN()-1, xer, fpe[1]/(corr*mpv_cosmics));
  }//if graph
  return {fp[1],fpe[1]};
}
////////////////////////////////////////////////////////////////////////////////

void fill_2d_graph(TGraph2D* graph, double x, double y, double z){
  graph->SetPoint(graph->GetN(), x, y, z/mpv_cosmics);
}

void fill_2d_hist(TH2D* h, double x, double y, double z){
  float Z = -1.;
  Z = z/mpv_cosmics;
  int bin = h->FindBin(x,y);
  float content = h->GetBinContent(bin);
  if( Z > content){
    h->SetBinContent(bin,Z);
  }
  return;
}

vector<double> ReadFit(TH1D* hdQds, TGraphErrors* graph, float x, float xer){
  TF1 *function;
  TList* funcs = hdQds->GetListOfFunctions();
  if(funcs->GetSize() == 0){
    #if verbose
    cout << "Bad fit for " << hdQds->GetName() << ". Ignoring." << endl;
    #endif
    return {-1,-1};
  }
  function = (TF1*)funcs->At(0);
  short color = function->GetLineColor();
  if(color != 600){
    #if verbose
    cout << "Bad fit for " << hdQds->GetName() << ". Ignoring." << endl;
    #endif
    return {-1,-1};
  }
  #if verbose
  cout << "Read fit: " << hdQds->GetName() << endl;
  #endif
  double MPV = function->GetMaximumX();
//  double MPV = function->GetParameter(1);
  double fpe = function->GetParError(1);
  if( graph != 0 ){
    //insert points in graph
    //insert points in graph
    graph->SetPoint(graph->GetN(), x, MPV/(corr*mpv_cosmics));
    graph->SetPointError(graph->GetN()-1, xer, fpe/(corr*mpv_cosmics));
  }//if graph
  return {MPV,fpe};
}


bool init_graph_gain(vector<TGraphErrors*> &mpv_field, map<int, vector<TGraphErrors*> > &mpv_field_ByLEMs, vector<TMultiGraph*> &mpv_field_AllLEMs, vector<int> &scan_nums_for_AllLEMs, string scan_type, int scan_num){

  mpv_field_AllLEMs.push_back(new TMultiGraph());mpv_field_AllLEMs.push_back(new TMultiGraph());mpv_field_AllLEMs.push_back(new TMultiGraph());
  mpv_field_AllLEMs[0]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_AllLEMs_view0").data());
  mpv_field_AllLEMs[1]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_AllLEMs_view1").data());
  mpv_field_AllLEMs[2]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_AllLEMs_summedviews").data());
      
  mpv_field.push_back(new TGraphErrors()); mpv_field.push_back(new TGraphErrors()); mpv_field.push_back(new TGraphErrors());
  mpv_field[0]->SetMarkerStyle(31);
  mpv_field[1]->SetMarkerStyle(31);
  mpv_field[2]->SetMarkerStyle(31);
  mpv_field[0]->SetMarkerSize(2);
  mpv_field[1]->SetMarkerSize(2);
  mpv_field[2]->SetMarkerSize(2);
  mpv_field[0]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_view0").data());
  mpv_field[1]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_view1").data());
  mpv_field[2]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_summed_views").data());
      
  for( auto lem : lems ){
    scan_nums_for_AllLEMs.push_back(scan_num);
    mpv_field_ByLEMs[lem] = {};
    mpv_field_ByLEMs[lem].push_back(new TGraphErrors()); mpv_field_ByLEMs[lem].push_back(new TGraphErrors()); mpv_field_ByLEMs[lem].push_back(new TGraphErrors());
    mpv_field_ByLEMs[lem][0]->SetMarkerStyle(31);
    mpv_field_ByLEMs[lem][1]->SetMarkerStyle(31);
    mpv_field_ByLEMs[lem][2]->SetMarkerStyle(31);
    mpv_field_ByLEMs[lem][0]->SetMarkerSize(2);
    mpv_field_ByLEMs[lem][1]->SetMarkerSize(2);
    mpv_field_ByLEMs[lem][2]->SetMarkerSize(2);
    mpv_field_ByLEMs[lem][0]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_view0").data());
    mpv_field_ByLEMs[lem][1]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_view1").data());
    mpv_field_ByLEMs[lem][2]->SetName(string("gain_"+scan_type+"_"+to_string(scan_num)+"_LEM_"+to_string(lem)+"_summed_views").data());
  }
  return true;
}

bool init_graph_purity_and_charging_up(vector<TGraphErrors> &graph,  map<int, vector<TGraphErrors> > &graph_ByLEMs, string name){
  
  if( dx > 100. ){
    #if verbose
    cout << "ERROR: can not have dx superior to total height of detector" << endl;
    #endif
    return false;
  }

  //multigraph for mpv, containing view 0 and view 1
  graph.push_back(TGraphErrors()); graph.push_back(TGraphErrors()); graph.push_back(TGraphErrors());
  graph[0].SetMarkerStyle(31);
  graph[1].SetMarkerStyle(31);
  graph[2].SetMarkerStyle(31);
  graph[0].SetMarkerStyle(2);
  graph[1].SetMarkerStyle(2);
  graph[2].SetMarkerStyle(2);
  graph[0].SetMarkerColor(kRed);
  graph[1].SetMarkerColor(kBlue);
  graph[0].SetName(string(name+"_view0").data());
  graph[1].SetName(string(name+"_view1").data());
  graph[2].SetName(string(name+"_summed_view").data());
  
  for( auto lem : lems ){
    graph_ByLEMs[lem] = {};
    graph_ByLEMs[lem].push_back(TGraphErrors()); graph_ByLEMs[lem].push_back(TGraphErrors()); graph_ByLEMs[lem].push_back(TGraphErrors());
    graph_ByLEMs[lem][0].SetMarkerStyle(31);
    graph_ByLEMs[lem][1].SetMarkerStyle(31);
    graph_ByLEMs[lem][2].SetMarkerStyle(31);
    graph_ByLEMs[lem][0].SetMarkerStyle(2);
    graph_ByLEMs[lem][1].SetMarkerStyle(2);
    graph_ByLEMs[lem][2].SetMarkerStyle(2);
    graph_ByLEMs[lem][0].SetMarkerColor(kRed);
    graph_ByLEMs[lem][1].SetMarkerColor(kBlue);
    graph_ByLEMs[lem][0].SetName(string(name+"_LEM_"+to_string(lem)+"_view0").data());
    graph_ByLEMs[lem][1].SetName(string(name+"_LEM_"+to_string(lem)+"_view1").data());
    graph_ByLEMs[lem][2].SetName(string(name+"_LEM_"+to_string(lem)+"_summed_view").data());
  }
  return true;
}

void fill_multigraph(map<int,vector<TGraphErrors*> > graphs_ByLems, vector<TMultiGraph*> &graph){

  if(MyColorPalette.size() < graphs_ByLems.size()){
    cout << "ERROR: not enoug hcolors in color palette" << endl;
    return;
  }
  int dcolor = (int)MyColorPalette.size()/graphs_ByLems.size();

  int i = 0;
  for(auto gr : graphs_ByLems){
    gr.second[0]->SetLineColor(MyColorPalette[i*dcolor]);
    gr.second[1]->SetLineColor(MyColorPalette[i*dcolor]);
    gr.second[2]->SetLineColor(MyColorPalette[i*dcolor]);
    gr.second[0]->SetMarkerColor(MyColorPalette[i*dcolor]);
    gr.second[1]->SetMarkerColor(MyColorPalette[i*dcolor]);
    gr.second[2]->SetMarkerColor(MyColorPalette[i*dcolor]);
    graph[0]->Add(new TGraphErrors(*gr.second[0]),"*");
    ((TGraphErrors*)graph[0]->GetListOfGraphs()->At(i))->SetTitle(string("LEM "+ to_string(gr.first)).data());
    graph[1]->Add(new TGraphErrors(*gr.second[1]),"*");
    ((TGraphErrors*)graph[1]->GetListOfGraphs()->At(i))->SetTitle(string("LEM "+ to_string(gr.first)).data());
    graph[2]->Add(new TGraphErrors(*gr.second[2]),"*");
    ((TGraphErrors*)graph[2]->GetListOfGraphs()->At(i))->SetTitle(string("LEM "+ to_string(gr.first)).data());
    i++;
  }
  return;
}

void sum_views(vector<TGraphErrors> &graph, map<int, vector<TGraphErrors> > &graph_ByLEMs){

  sum_views_singlegraph(graph);
  
  for( auto lem : lems){
    if( graph_ByLEMs[lem][0].GetN() > 2 and graph_ByLEMs[lem][1].GetN() > 2 ){sum_views_singlegraph(graph_ByLEMs[lem]);}
  }//for lem
  return;
}

void sum_views(vector<TGraphErrors*> &graph, map<int, vector<TGraphErrors*> > &graph_ByLEMs){

  sum_views_singlegraph(graph);
  
  for( auto lem : lems ){
    if( graph_ByLEMs[lem][0]->GetN() > 2 and graph_ByLEMs[lem][1]->GetN() > 2 ){sum_views_singlegraph(graph_ByLEMs[lem]);}
  }//for lem
  return;
}

void sum_views_singlegraph(vector<TGraphErrors> &graph){
  int biggest_graph = 0;
  int smallest_graph = 1;
  if( graph[0].GetN() < graph[1].GetN() ){
    biggest_graph = 1;
    smallest_graph = 0;
  }
  for( int i = 0; i < graph[biggest_graph].GetN(); i++ ){
    double x1 = 0;
    double y1 = 0;
    double ex1 = 0;
    double y3 = 0;
    double ey3 = 0;
    double x2 = 0;
    double y2 = 0;
    graph[biggest_graph].GetPoint(i,x1,y1);
    double ey1 = graph[biggest_graph].GetErrorY(i);
    for( int j = 0; j < graph[biggest_graph].GetN(); j++ ){
      graph[smallest_graph].GetPoint(j,x2,y2);
      double ey2 = graph[smallest_graph].GetErrorY(j);
      if( x2 == x1 ){
        y3 = y1+y2;
        ey3 = TMath::Sqrt(ey1*ey1 + ey2*ey2);
        break;
      }
    }
    if( y3 != y1 and y3 != 0 ){
      graph[2].SetPoint(graph[2].GetN(),x1,y3);
      graph[2].SetPointError(graph[2].GetN()-1,0,ey3);
    }
  }
}

void sum_views_singlegraph(vector<TGraphErrors*> &graph){
  int biggest_graph = 0;
  int smallest_graph = 1;
  if( graph[0]->GetN() < graph[1]->GetN() ){
    biggest_graph = 1;
    smallest_graph = 0;
  }
  for( int i = 0; i < graph[biggest_graph]->GetN(); i++ ){
    double x1 = 0;
    double y1 = 0;
    double ex1 = 0;
    double y3 = 0;
    double ey3 = 0;
    double x2 = 0;
    double y2 = 0;
    graph[biggest_graph]->GetPoint(i,x1,y1);
    double ey1 = graph[biggest_graph]->GetErrorY(i);
    for( int j = 0; j < graph[biggest_graph]->GetN(); j++ ){
      graph[smallest_graph]->GetPoint(j,x2,y2);
      double ey2 = graph[smallest_graph]->GetErrorY(j);
      if( x2 == x1 ){
        y3 = y1+y2;
        ey3 = TMath::Sqrt(ey1*ey1 + ey2*ey2);
        break;
      }
    }
    if( y3 != y1 and y3 != 0 ){
      graph[2]->SetPoint(graph[2]->GetN(),x1,y3);
      graph[2]->SetPointError(graph[2]->GetN()-1,0,ey3);
    }
  }
}

double RFromBirk(double drift, double MPV, bool convert){
  double A = 0.8;
  double k = 0.0486; //kV/cm
  double rho = 1.3954; //g/cm3
  if(MPV > 0){
    if(convert){
      MPV = 1./(1.6*A/(0.236*mpv_cosmics) - k/0.5);
    }
    return A/(1+MPV*k/(rho*drift));
  }
  else if(mpv_cosmics < 0){
    if(!load_cosmics()){cout << "ERROR in RFromBirk: can't load cosmics" << endl; return -1.;}
  }
  double E = 1./(1.6*A/(0.236*mpv_cosmics) - k/0.5);
  return A/(1+E*k/drift);
}

double MPVs_from_theory(double ds, double drift){

  double mpv;
  
  double a = 0.19559; //no units
  double k = 3.; //no units
  double x0 = 0.2; //no units
  double x1 = 3.; //no units
  double d0 = 0.; //no units
  double C = 5.2146; //no units
  double M = 105.658;//MeV/c2
  double I = 1.88e-4; //MeV
  double A = 39.948; //no units
  double Z = 18.; //no units
  double q = 1.; //no units
  double K = 0.307075; //MeV/mol cm2
  double me = 0.511; //MeV/c**2
  double rho = 1.3954; //g/cm3
  double E = 4000.;//4 GeV;/c**2
  double p = TMath::Sqrt(E*E-M*M); 
  double beta = p/E;
  double gamma = 1/TMath::Sqrt(1-beta*beta);
  double x = TMath::Log10(beta*gamma);
  double dbg = 0;
  
  if(x > x1){dbg = 2*TMath::Log(10)*x-C;}
  else if(x < x1 and x > x0){dbg = 2*TMath::Log(10)*x-C + a*TMath::Power((x1-x),k);}
  else if(x < x0){dbg = d0*TMath::Power(10,(2*x-x0));}
  double eta = 0.5*K*(Z/A)*(rho*ds)/(beta*beta);
  double MPV = (K*Z*rho/(2*A*beta*beta)) * (TMath::Log(2*me*rho*K*Z*gamma*gamma/(2*A*I*I)) + 0.2-beta*beta-dbg) + K*Z*rho*TMath::Log(ds)/(2*A*beta*beta);
  
  double R =  RFromBirk(drift, MPV); //all in cm, MeV, g and kV
  
  return 1.6*MPV*R/0.236;
}

double e_lifetime(TGraphErrors *graph){
  double G0=0.;
  double x0=0;
  double G=100.;
  double x=0;
  for(int i = 0; i < graph->GetN(); i++){
    double dummyx, dummyG;
    graph->GetPoint(i, dummyx, dummyG);
    if(dummyG > G0){
      G0 = dummyG;
      x0 = dummyx;
    }
    if(dummyG < G){
      G = dummyG;
      x = dummyx;
    }
  }
  TF1 *myfit = new TF1("myfit","[0]*exp(-(x-[2])/(0.16*[1]))",x0,x);
  myfit->SetLineColor(graph->GetMarkerColor());
  myfit->SetParameter(0,G0);
  myfit->SetParameter(1,5000.);
  myfit->SetParameter(2,x0);
  graph->Fit("myfit");
  //v elec in LAr at 500V/cm ~ 1.6 mm/micro sec
  return myfit->GetParameter(1);
}

double get_eff(double extraction, double amplification, double induction){
  double eff = -1.;
  
  if(amplification < 0 and induction < 0){
    if(extraction > 100){extraction = extraction/1000.;}
    if(h_GushinEff.GetEntries() == 0){load_Gushin_Eff();}
    int bin = h_GushinEff.FindBin(extraction);
    eff = h_GushinEff.GetBinContent(bin);
    return eff;
  }
  
  //convert V/cm or kV/cm to V
  amplification = amplification*100.;
  extraction = extraction*1000.;
  induction = induction*200.;
  if (induction == 5){
    induction = 4.999;
  }
  
  
  if(h_ExtrEff_vs_LemExtr.GetEntries() == 0){
    if(!load_extr_eff_simu_graphs()){return -1;}
  }
  if(h_IndEff_vs_LemInd.GetEntries() == 0){
    if(!load_ind_eff_simu_graphs()){return -1;}
  }
  if(extraction > 0){
    int bin = h_ExtrEff_vs_LemExtr.FindBin(amplification,extraction);
    eff = h_ExtrEff_vs_LemExtr.GetBinContent(bin);
  }
  else{
    eff = 1.;
  }
  if(induction > 0){
    int bin = h_IndEff_vs_LemInd.FindBin(amplification,induction);
    eff = eff*h_IndEff_vs_LemInd.GetBinContent(bin);
  }
  return eff;
}

pair<TH2D*, TGraph*> normalise_gain_histo_2d(TH2D* h, string scan_type){
  TGraph *ratio = NULL;
  TH2D *corresponding_eff = NULL;
  if(scan_type != "Amplification_Induction" and scan_type != "Amplification_Extraction"){return make_pair(corresponding_eff, ratio);}
  vector<int> run_list;
  vector<pair<float,float> > field;
  vector<string> const_fields_names;
  vector<float> const_fields_values;
  if(!load_runs_2D(scan_type, atoi(h->GetTitle()), run_list, field, const_fields_names, const_fields_values)){return make_pair(corresponding_eff, ratio);}
  ratio = new TGraph();
  ratio->SetName(string(string(h->GetName())+"_ratio").data());
  ratio->SetTitle(string(string(ratio->GetName()) + ";bin;normalized gain/normalized simu").data());
  string field1 = scan_type.substr(0,scan_type.find_first_of("_"));
  string field2 = scan_type.substr(scan_type.find_first_of("_")+1);
  corresponding_eff = (TH2D*)h->Clone();
  corresponding_eff->SetName(string(string(h->GetName())+"_corresponding_eff").data());
  double amplifield = -1.;
  double inducorextrfield = -1.;
  double gain = -1.;
  double maxeff = -1.;
  double maxgain = -1.;
  double eff = -1.;
  int bin = -2.;
  double norm = -1.;
  int ratio_x = 0;
  for(int ix = corresponding_eff->GetNbinsX()-1; ix >= 0; ix--){
    amplifield = corresponding_eff->GetXaxis()->GetBinCenter(ix) - (corresponding_eff->GetXaxis()->GetBinCenter(2)-corresponding_eff->GetXaxis()->GetBinCenter(1))/2.;
    for(int iy = corresponding_eff->GetNbinsY()-1 ; iy >= 0; iy--){
      inducorextrfield = corresponding_eff->GetYaxis()->GetBinCenter(iy) - (corresponding_eff->GetYaxis()->GetBinCenter(2)-corresponding_eff->GetYaxis()->GetBinCenter(1))/2.;
      bin = corresponding_eff->GetBin(ix,iy);
      gain = h->GetBinContent(bin);
      if(gain < 1e-5){continue;}
      if(scan_type == "Amplification_Induction"){eff = get_eff(-1.,amplifield,inducorextrfield);}
      else if(scan_type == "Amplification_Extraction"){eff = get_eff(inducorextrfield,amplifield,-1.);}
      if(maxgain < gain){maxgain = gain;}
      if(maxeff < eff){maxeff = eff;}
      corresponding_eff->SetBinContent(bin,eff);
    }
  }
  for(int ix = corresponding_eff->GetNbinsX()-1; ix >= 0; ix--){
    amplifield = corresponding_eff->GetXaxis()->GetBinCenter(ix) - (corresponding_eff->GetXaxis()->GetBinCenter(2)-corresponding_eff->GetXaxis()->GetBinCenter(1))/2.;
    for(int iy = corresponding_eff->GetNbinsY()-1 ; iy >= 0; iy--){
      inducorextrfield = corresponding_eff->GetYaxis()->GetBinCenter(iy) - (corresponding_eff->GetYaxis()->GetBinCenter(2)-corresponding_eff->GetYaxis()->GetBinCenter(1))/2.;
      bin = corresponding_eff->GetBin(ix,iy);
      gain = h->GetBinContent(bin)/maxgain;
      eff =  corresponding_eff->GetBinContent(bin)/maxeff; 
      corresponding_eff->SetBinContent(bin,eff);
      if(gain > 0){
        if(scan_type == "Amplification_Extraction"){
          if(inducorextrfield >= 1.6){
            ratio->SetPoint(ratio->GetN(),ratio_x+1, gain/eff);
            ratio_x++;
          }
        }
        else{
          ratio->SetPoint(ratio->GetN(),ratio_x+1, gain/eff);
          ratio_x++;
        }
      }
    }
  }
  return make_pair(corresponding_eff, ratio);
}

void get_eff_multigraphs(TMultiGraph *mg, string scan_type, vector<int> scan_nums, TFile *ofile, string eff_to_use, string draw){
  double max = 0;
  if(scan_type == "Amplification"){
    #if verbose
    cout << "get_eff_multigraphs: No normalization implemented for amplification scan yet" << endl;
    #endif
  }
  if(MyColorPalette.size() < (unsigned int)mg->GetListOfGraphs()->GetSize()){
    cout << "ERROR: not enough hcolors in color palette" << endl;
    return;
  }
  int dcolor = (int)MyColorPalette.size()/mg->GetListOfGraphs()->GetSize();
  vector<TGraphErrors*> gr_norm;
  vector<TGraphErrors*> mg_eff;
  vector<TGraphErrors*> mg_ratio;
  vector<TH1D*> vec_distri_above_16;
  vector<TH1D*> vec_distri;
  TCanvas *can_ratio = new TCanvas(string("can_" + string(mg->GetName()) + "_ratio").data(),string("can_" + string(mg->GetName()) + "_ratio").data(),1500,1500/golden_ratio);
  TCanvas *can_distri = new TCanvas(string("can_" + string(mg->GetName()) + "_distri").data(),string("can_" + string(mg->GetName()) + "_distri").data(),1500,1500/golden_ratio);
  TCanvas *can_distri_16 = new TCanvas(string("can_" + string(mg->GetName()) + "_distri_16").data(),string("can_" + string(mg->GetName()) + "_distri_16").data(),1500,1500/golden_ratio);
  if(eff_to_use == "Gushin"){
    can_ratio->SetName(string(string(can_ratio->GetName())+"_Gushin").data());
    can_distri->SetName(string(string(can_distri->GetName())+"_Gushin").data());
    can_distri_16->SetName(string(string(can_distri_16->GetName())+"_Gushin").data());
  }
  TCanvas *can_both = new TCanvas(string("can_" + string(mg->GetName()) + "_both").data(),string("can_" + string(mg->GetName()) + "_both").data(),1500,1500/golden_ratio);
  
  bool first = true;
  
  for(int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++){
    gr_norm.push_back((TGraphErrors*)mg->GetListOfGraphs()->At(i));
    mg_eff.push_back((TGraphErrors*)mg->GetListOfGraphs()->At(i));
    mg_ratio.push_back((TGraphErrors*)mg->GetListOfGraphs()->At(i));
    vec_distri_above_16.push_back(new TH1D(string(string(mg_eff.back()->GetName())+"_distri_16").data(),string(mg_eff.back()->GetTitle()).data(),50,0,2));
    vec_distri.push_back(new TH1D(string(string(mg_eff.back()->GetName())+"_distri").data(),string(mg_eff.back()->GetTitle()).data(),50,0,2));
    gr_norm.back()->SetName(string(string(gr_norm.back()->GetName())+"_norm").data());
    mg_eff.back()->SetName(string(string(mg_eff.back()->GetName())+"_eff").data());
    mg_ratio.back()->SetName(string(string(mg_ratio.back()->GetName())+"_ratio").data());
    if(eff_to_use == "Gushin"){
      mg_eff.back()->SetName(string(string(mg_eff.back()->GetName())+"_Gushin").data());
      mg_ratio.back()->SetName(string(string(mg_ratio.back()->GetName())+"_Gushin").data());
    }
  }
    
  for(unsigned int i = 0; i<mg_eff.size(); i++){
    gr_norm[i]->SetMarkerColor(MyColorPalette[i*dcolor]);
    mg_eff[i]->SetMarkerColor(MyColorPalette[i*dcolor]);
    mg_ratio[i]->SetMarkerColor(MyColorPalette[i*dcolor]);
    gr_norm[i]->SetLineColor(MyColorPalette[i*dcolor]);
    mg_eff[i]->SetLineColor(MyColorPalette[i*dcolor]);
    mg_ratio[i]->SetLineColor(MyColorPalette[i*dcolor]);
    vec_distri_above_16[i]->SetLineColor(MyColorPalette[i*dcolor]);
    vec_distri[i]->SetLineColor(MyColorPalette[i*dcolor]);
    gr_norm[i]->SetMarkerStyle(22);
    mg_eff[i]->SetMarkerStyle(29);
    mg_ratio[i]->SetMarkerStyle(29);
    
    vector<int> run_list;
    vector<float> dummy;
    vector<string> const_fields_names;
    vector<float> const_fields_values;
    if(!load_runs(scan_type, scan_nums[i], run_list, dummy, const_fields_names, const_fields_values)){return;}
    double amplifield = -1.;
    double extracfield = -1.;
    double inducfield = -1.;
    double driftfield = -1.;
    
    for( unsigned int j = 0; j < const_fields_values.size(); j++){
      if( const_fields_names.at(j) == "Extraction" ){extracfield = const_fields_values.at(j);}
      else if( const_fields_names.at(j) == "Induction" ){inducfield = const_fields_values.at(j);}
      else if( const_fields_names.at(j) == "Amplification" ){amplifield = const_fields_values.at(j);}
      else if( const_fields_names.at(j) == "Drift" ){driftfield = const_fields_values.at(j);}
    }
    
    double gain = 0;
    double field = 0;
    double norm = 0;
    double eff = 0;
    double max_eff = 0;
    double max_gain = 0;
    for(int j = 0; j < mg_eff[i]->GetN(); j++){
      mg_eff[i]->GetPoint(j,field,gain);
      if(scan_type == "Drift"){eff = RFromBirk(field, gain*mpv_cosmics, true);}
      if(scan_type == "Induction"){eff = get_eff(-1,amplifield,field);}
      if(scan_type == "Extraction"){
        if(eff_to_use == "simu"){
          eff = get_eff(field,amplifield,-1);
        }
        else if(eff_to_use == "Gushin"){
          eff = get_eff(field);
        }
      }
      if(max_eff < eff){max_eff = eff;}
      if(max_gain < gain){max_gain = gain;}
      mg_eff[i]->SetPoint(j,field,eff);
    }
    for(int j = 0; j < mg_eff[i]->GetN(); j++){
      gr_norm[i]->GetPoint(j,field,gain);
      gr_norm[i]->SetPoint(j,field,gain/max_gain);
      mg_eff[i]->GetPoint(j,field,eff);
      mg_eff[i]->SetPoint(j,field,eff/max_eff);
      mg_ratio[i]->SetPoint(j,field,gain*max_eff/(eff*max_gain));
      vec_distri[i]->Fill(gain*max_eff/(eff*max_gain));
      if(field >= 1.6){
        vec_distri_above_16[i]->Fill(gain*max_eff/(eff*max_gain));
      }
    }
    if(first){
      can_ratio->cd();
      mg_ratio[i]->Draw("AP");
      mg_ratio[i]->SetMinimum(0);
      mg_ratio[i]->Draw("AP");
      can_ratio->Update();
      
      can_distri->cd();
      vec_distri[i]->SetMaximum(20);
      vec_distri[i]->Draw();
      
      can_distri_16->cd();
      vec_distri_above_16[i]->SetMaximum(20);
      vec_distri_above_16[i]->Draw();
      
      can_both->cd();
      gr_norm[i]->Draw("AP");
      gr_norm[i]->SetMinimum(0);
      gr_norm[i]->Draw("AP");
      can_both->Update();
      mg_eff[i]->Draw("SAME P");
      
      first = false;
    }
    else{
      can_ratio->cd();
      mg_ratio[i]->Draw("SAME P");
      
      can_distri->cd();
      vec_distri[i]->Draw("SAME");
      
      can_distri_16->cd();
      vec_distri_above_16[i]->Draw("SAME");
      
      can_both->cd();
      gr_norm[i]->Draw("SAME P");
      mg_eff[i]->Draw("SAME P");
    }
  }
  
  ofile->cd();
  TLegend *leg_ratio = can_ratio->BuildLegend(.7,.1,.9,.3);
  can_ratio->Write();
  if(draw != ""){
    can_ratio->SaveAs(string(draw + string(can_ratio->GetName())+".png").data());
  }
  delete can_ratio;
  
  ofile->cd();
  TLegend *leg_distri = can_distri->BuildLegend(.7,.6,.9,.9);
  can_distri->Write();
  if(draw != ""){
    can_distri->SaveAs(string(draw + string(can_distri->GetName())+".png").data());
  }
  delete can_distri;
  
  ofile->cd();
  TLegend *leg_distri_16 = can_distri_16->BuildLegend(.7,.6,.9,.9);
  can_distri_16->Write();
  if(draw != ""){
    can_distri_16->SaveAs(string(draw + string(can_distri_16->GetName())+".png").data());
  }
  delete can_distri_16;
  
  ofile->cd();
  TLegend *leg_both = can_both->BuildLegend(.7,.6,.9,.9);
  can_both->Write();
  if(draw != ""){
    can_both->SaveAs(string(draw + string(can_both->GetName())+".png").data());
  }
  delete can_both;
  
  gr_norm.clear();
  vec_distri.clear();
  vec_distri_above_16.clear();
  mg_eff.clear();
  mg_ratio.clear();
  return;
}

void get_eff_graphs(TGraphErrors *gr, string scan_type, int scan_num, TFile *ofile, string eff_to_use, string draw){
  if(scan_type == "Amplification"){return;}
  
  TGraphErrors* gr_eff = (TGraphErrors*)gr->Clone();
  TGraphErrors* gr_ratio = (TGraphErrors*)gr->Clone();
  gr_eff->SetName(string(string(gr->GetName())+"_effs").data());
  gr_ratio->SetName(string(string(gr->GetName())+"_ratio").data());
  TH1D *distri_above_16 = new TH1D(string(string(gr->GetName())+"_distri_16").data(),"",50,0,2);
  TH1D *distri = new TH1D(string(string(gr->GetName())+"_distri").data(),"",50,0,2);

  if(eff_to_use == "Gushin"){
    gr->SetName(string(string(gr->GetName())+"_Gushin").data());
    gr_eff->SetName(string(string(gr_eff->GetName())+"_Gushin").data());
    gr_ratio->SetName(string(string(gr_ratio->GetName())+"_Gushin").data());
    distri_above_16->SetName(string(string(distri_above_16->GetName())+"_Gushin").data());
    distri->SetName(string(string(distri->GetName())+"_Gushin").data());
  }
  vector<int> run_list;
  vector<float> dummy;
  vector<string> const_fields_names;
  vector<float> const_fields_values;
  if(!load_runs(scan_type, scan_num, run_list, dummy, const_fields_names, const_fields_values)){return;}
  double amplifield = -1.;
  double extracfield = -1.;
  double inducfield = -1.;
  double driftfield = -1.;
  
  for( unsigned int j = 0; j < const_fields_values.size(); j++){
    if( const_fields_names.at(j) == "Extraction" ){extracfield = const_fields_values.at(j);}
    else if( const_fields_names.at(j) == "Induction" ){inducfield = const_fields_values.at(j);}
    else if( const_fields_names.at(j) == "Amplification" ){amplifield = const_fields_values.at(j);}
    else if( const_fields_names.at(j) == "Drift" ){driftfield = const_fields_values.at(j);}
  }
  
  double gain = 0;
  double field = 0;
  double norm = 0;
  double eff = 0;
  double max_eff = 0;
  double max_gain = 0;
  for(int j = 0; j < gr_eff->GetN(); j++){
    gr_eff->GetPoint(j,field,gain);
    if(scan_type == "Drift"){eff = RFromBirk(field, gain*mpv_cosmics, true);}
    if(scan_type == "Induction"){eff = get_eff(-1,amplifield,field);}
    if(scan_type == "Extraction"){
      if(eff_to_use == "simu"){
        eff = get_eff(field,amplifield,-1);
      }
      else if(eff_to_use == "Gushin"){
        eff = get_eff(field);
      }
    }
    if(max_eff < eff){max_eff = eff;}
    if(max_gain < gain){max_gain = gain;}
    gr_eff->SetPoint(j,field,eff);
  }
  for(int j = 0; j < gr_eff->GetN(); j++){
    gr->GetPoint(j,field,gain);
    gr->SetPoint(j,field,gain/max_gain);
    gr_eff->GetPoint(j,field,eff);
    gr_eff->SetPoint(j,field,eff/max_eff);
    gr_ratio->SetPoint(j,field,gain*max_eff/(eff*max_gain));
    distri->Fill(gain*max_eff/(eff*max_gain));
    if(field > 1.8){
      distri_above_16->Fill(gain*max_eff/(eff*max_gain));
    }
  }
  
  if(gr_eff == 0){return;}
  gr_ratio->SetMarkerColor(kRed);
  gr_ratio->SetMarkerStyle(29);
  gr_ratio->SetLineColor(kRed);
  gr_ratio->Draw("AP*");
  gr_ratio->SetMinimum(0);
  gr_ratio->SetTitle("");
  gr_ratio->Draw("AP*");
  gPad->Update();
  gPad->SetName(string(string(gr_ratio->GetName())+"_pad").data());
  ofile->cd();
  gPad->Write();
  if(draw != ""){
    gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
  }
  delete gPad;
  
  gr->SetMarkerColor(kGreen);
  gr->SetLineColor(kGreen);
  gr->SetMarkerStyle(29);
  gr->Draw("AP*");
  gr->SetMinimum(0);
  gr->SetTitle("");
  gr->Draw("AP*");
  gr_eff->SetMarkerColor(kRed);
  gr_eff->SetLineColor(kRed);
  gr_eff->Draw("P* SAME");
  gPad->Update();
  gPad->SetName(string(string(gr->GetName())+"_both_pad").data());
  ofile->cd();
  gPad->Write();
  if(draw != ""){
    gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
  }
  delete gPad;
  
  distri->Draw();
  gPad->SetName(string(string(distri->GetName())+"_pad").data());
  ofile->cd();
  gPad->Write();
  if(draw != ""){
    gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
  }
  delete gPad;
  
  distri_above_16->Draw();
  gPad->SetName(string(string(distri_above_16->GetName())+"_pad").data());
  ofile->cd();
  gPad->Write();
  if(draw != ""){
    gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
  }
  delete gPad;
  delete distri;
  delete distri_above_16;
  delete gr_eff;
  
  return;
}

TMultiGraph* normalise_gain_graph(TMultiGraph *mg, string scan_type, vector<string> effs_to_use){
  TMultiGraph* chien = NULL;
  TMultiGraph *mg_normalized = (TMultiGraph*)mg->Clone();
  mg_normalized->SetName(string(string(mg->GetName())+"_norm").data());
  for(int i = 0; i < mg_normalized->GetListOfGraphs()->GetSize(); i++){
    TGraphErrors* gr = (TGraphErrors*)mg->GetListOfGraphs()->At(i);
    TGraphErrors* gr_normalized = (TGraphErrors*)mg_normalized->GetListOfGraphs()->At(i);
    gr_normalized->SetName(string(string(gr->GetName())+"_effs").data());
    vector<int> run_list;
    vector<float> field;
    vector<string> const_fields_names;
    vector<float> const_fields_values;
    if(!load_runs(scan_type, atoi(gr->GetTitle()), run_list, field, const_fields_names, const_fields_values)){continue;}
    double amplifield = -1.;
    double extracfield = -1.;
    double inducfield = -1.;
    double driftfield = -1.;
    double gain = -1.;
    double eff = 1.;
    
    for( unsigned int j = 0; j < const_fields_values.size(); j++){
      if( const_fields_names.at(j) == "Extraction" ){
        extracfield = const_fields_values.at(j);
      }
      else if( const_fields_names.at(j) == "Induction" ){
        inducfield = const_fields_values.at(j);
      }
      else if( const_fields_names.at(j) == "Amplification" ){
        amplifield = const_fields_values.at(j);
      }
      else if( const_fields_names.at(j) == "Drift" ){
        driftfield = const_fields_values.at(j);
      }
    }
  
    for(int i=0; i<gr_normalized->GetN(); i++){
      if(scan_type == "Amplification" or scan_type == "All"){
        gr_normalized->GetPoint(i,amplifield,gain);
      }
      else if(scan_type == "Induction"){
        gr_normalized->GetPoint(i,inducfield,gain);
      }
      else if(scan_type == "Extraction"){
        gr_normalized->GetPoint(i,extracfield,gain);
      }
      else if(scan_type == "Drift" or scan_type == "All"){
        gr_normalized->GetPoint(i,driftfield,gain);
      }
      for(auto eff_to_use : effs_to_use){
        if(eff_to_use == "Gushin"){eff = eff*get_eff(extracfield);}
        else if(eff_to_use == "Extr"){eff = eff*get_eff(extracfield,amplifield,-1.);}
        else if(eff_to_use == "Ind"){eff = eff*get_eff(-1.,amplifield,inducfield);}
        else if(eff_to_use == "Drift"){eff = eff*RFromBirk(driftfield, gain*mpv_cosmics, true);}
        else{cout << "WARNING in normalise_gain_graph: unknown eff " << eff_to_use << endl;}
        gr_normalized->SetName(string(string(gr_normalized->GetName())+"_"+eff_to_use).data());
      }
      if (eff <= 0 or eff > .999999){
        if(scan_type == "Amplification" or scan_type == "All"){
          gr_normalized->SetPoint(i,amplifield,-1);
        }
        else if(scan_type == "Induction"){
          gr_normalized->SetPoint(i,inducfield,-1);
        }
        else if(scan_type == "Extraction"){
          gr_normalized->SetPoint(i,extracfield,-1);
        }
        else if(scan_type == "Drift"){
          gr_normalized->SetPoint(i,driftfield,-1);
        }
      }
      else{
        if(scan_type == "Amplification" or scan_type == "All"){
          gr_normalized->SetPoint(i,amplifield,gain/eff);
        }
        else if(scan_type == "Induction"){
          gr_normalized->SetPoint(i,inducfield,gain/eff);
        }
        else if(scan_type == "Extraction"){
          gr_normalized->SetPoint(i,extracfield,gain/eff);
        }
        else if(scan_type == "Drift"){
          gr_normalized->SetPoint(i,driftfield,gain/eff);
        }
      }
      eff = 1.;
    }
  }
  return mg_normalized;
}

double GetMaxOFMapMapTH1D(map<double,map<int,TH1D> > MapMapTH1D){
  double max = -1;
  for(auto MapTH1D : MapMapTH1D){
    for(auto h : MapTH1D.second){
      double tmpmax = h.second.GetMaximum();
      if(tmpmax > max){max = tmpmax;}
    }
  }
  return max;
}

void draw_gain_dQds(map<double,map<int,TH1D> > mpv_field, string scan_type, int scan_num, TFile *ofile, string path, int ilem){

  if(MyColorPalette.size() < mpv_field.size()){
    cout << "ERROR: not enough hcolors in color palette" << endl;
    return;
  }
  int dcolor = (int)MyColorPalette.size()/mpv_field.size();
  string name = scan_type + "_" + to_string(scan_num);
  if(ilem != -1){
    name = name + "_" + to_string(ilem);
  }
  TLegend *myleg = new TLegend(.8,.6,.9,.9);
  TCanvas *mycan = new TCanvas(("can_"+name).data(),name.data(),1500,1500/golden_ratio);
  mycan->SetGrid();
  
  double max = GetMaxOFMapMapTH1D(mpv_field);
  bool first = true;
  int color = 0;
  for(auto field : mpv_field){
    int col;
    int first_run = true;
    for(auto run : field.second){
      mpv_field[field.first][run.first].SetLineWidth(2);
      mpv_field[field.first][run.first].SetLineColor(MyColorPalette[color]);
      mpv_field[field.first][run.first].SetMarkerColor(MyColorPalette[color]);
      mpv_field[field.first][run.first].SetMarkerStyle(kFullTriangleDown);
      mpv_field[field.first][run.first].SetName(to_string_with_precision(field.first,3).data());
      mpv_field[field.first][run.first].SetTitle(to_string_with_precision(field.first,3).data());
      mpv_field[field.first][run.first].SetMarkerStyle(kFullTriangleDown);
      if(first){
        mpv_field[field.first][run.first].GetXaxis()->SetNdivisions(50);
        mpv_field[field.first][run.first].GetXaxis()->SetLabelSize(.02);
        mpv_field[field.first][run.first].SetMaximum(max*1.1);
        mpv_field[field.first][run.first].SetMinimum(0);
        mpv_field[field.first][run.first].GetXaxis()->SetTitle("dQ/ds (f/cm)");
        mpv_field[field.first][run.first].GetYaxis()->SetTitle("# hits (normalised)");
        mpv_field[field.first][run.first].Draw();
        myleg->AddEntry(&mpv_field[field.first][run.first]);
        first = false;
        first_run = false;
      }
      else if(first_run){
        myleg->AddEntry(&mpv_field[field.first][run.first]);
        first_run = false;
      }
      mpv_field[field.first][run.first].Draw("SAME");
    }
    color+=dcolor;
  }
  ofile->cd();
  myleg->Draw();
  mycan->Write();
  if(path != ""){
    mycan->SaveAs(string(path+string(mycan->GetName())+".png").data());
  }
  delete mycan;
  return;
}

//void draw_gain_multigraph(TMultiGraph *mg, string scan_type, TFile *ofile, bool draw_leg, string draw){
//  if(!mg->GetListOfGraphs()){
//    #if verbose
//    cout << "draw_gain_multigraph: " << mg->GetName() << " is empty " << endl;
//    #endif
//    return;
//  }
//  mg->Draw("A pmc");
//  mg->SetMinimum(0);
//  gPad->Update();
//  ofile->cd();
//  mg->Write();
//  if(draw_leg){
//    TLegend *leg = gPad->BuildLegend(.65,.1,.9,.17);
//    leg->SetNColumns(mg->GetListOfGraphs()->GetSize());
//    leg->SetFillColorAlpha(0,0);
//    for( int i=0; i<leg->GetListOfPrimitives()->GetSize(); i++ ){
//      TLegendEntry *legentry = (TLegendEntry*)leg->GetListOfPrimitives()->At(i);
//      legentry->SetOption("p");
//    }
//    gPad->SetName(string(string(mg->GetName())+"_pad").data());
//    ofile->cd();
//    gPad->Write();
//    if(draw != ""){
//      gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
//    }
//    delete leg;
//  }
//  else{
//    gPad->SetName(string(string(mg->GetName())+"_pad").data());
//    ofile->cd();
//    gPad->Write();
//    if(draw != ""){
//      gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
//    }
//  }
//  delete gPad;
//  vector<pair<TMultiGraph*,double> > vec_eff_graphs_and_max = {};
//  vec_eff_graphs_and_max.push_back(get_eff_multigraphs(mg, scan_type));
//  if(scan_type == "Extraction"){
//    vec_eff_graphs_and_max.push_back(get_eff_multigraphs(mg, scan_type, "Gushin"));
//  }
//  for(auto eff_graphs_and_max : vec_eff_graphs_and_max){
//    if(!eff_graphs_and_max.first){continue;}
//    TMultiGraph* eff_graphs = eff_graphs_and_max.first;
//    double max = eff_graphs_and_max.second;
//    for(int i=0; i<eff_graphs->GetListOfGraphs()->GetSize(); i++){
//      ((TGraphErrors*)mg->GetListOfGraphs()->At(i))->Draw("AP");
//      ((TGraphErrors*)mg->GetListOfGraphs()->At(i))->SetMinimum(0);
//      ((TGraphErrors*)mg->GetListOfGraphs()->At(i))->SetMaximum(max*1.1);
//      ((TGraphErrors*)mg->GetListOfGraphs()->At(i))->SetTitle("");
//      ((TGraphErrors*)mg->GetListOfGraphs()->At(i))->Draw("AP");
//      gPad->Update();
//      ((TGraphErrors*)eff_graphs->GetListOfGraphs()->At(i))->SetMarkerColor(kRed);
//      ((TGraphErrors*)eff_graphs->GetListOfGraphs()->At(i))->SetMarkerStyle(29);
//      ((TGraphErrors*)eff_graphs->GetListOfGraphs()->At(i))->SetLineColor(kRed);
//      ((TGraphErrors*)eff_graphs->GetListOfGraphs()->At(i))->Draw("SAME P");
//      gPad->SetName(string(string(eff_graphs->GetName())+"_pad").data());
//      ofile->cd();
//      gPad->Write();
//      if(draw != ""){
//        gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
//      }
//    }
//    delete gPad;
//    delete eff_graphs;
//  }
//  delete mg;
//  return;
//}

void draw_gain_graph(TGraphErrors *gr, string scan_type, int scan_num, TFile *ofile, string draw){
//  if(scan_type == "Amplification" and scan_num == 100){
//    if(Gain_graph_3L.GetN() == 0){load_fit_3L();}
//    Gain_graph_3L.Draw("AP*");
//    Gain_graph_3L.SetMinimum(0);
//    Gain_graph_3L.Draw("AP*");
//    gPad->Update();
//    gPad->Modified();
//    ((TF1*)Gain_graph_3L.GetListOfFunctions()->At(0))->Draw("SAME");
//    gr->Draw("SAME P*");
//  }
//  else{
    gr->Draw("AP*");
    gr->SetMinimum(0);
    gr->Draw("AP*");
    gPad->Update();
    gPad->Modified();
//  }
  ofile->cd();
  gr->Write();
  gPad->SetName(string(string(gr->GetName())+"_pad").data());
  ofile->cd();
  gPad->Write();
  if(draw != ""){
    gPad->SaveAs(string(draw + string(gPad->GetName())+".png").data());
  }
  delete gPad;
  if(scan_type == "Amplification"){delete gr; return;}
  get_eff_graphs(gr, scan_type, scan_num, ofile, "simu", draw);
  if(scan_type == "Extraction"){
    get_eff_graphs(gr, scan_type, scan_num, ofile, "Gushin", draw);
  }
  delete gr;
  return;
}

void draw_gain_graph_superposed_lems(TMultiGraph* gr, string scan_type, vector<int> scan_nums, TFile *ofile, string draw){

  TCanvas *mycan = new TCanvas(("can_"+string(gr->GetName())).data(),string(gr->GetName()).data(),1500,1500/golden_ratio);
  mycan->SetGrid();
  
  double max = 1.1;
  bool first = true;
  int color = 0;
  int col;
  int first_run = true;
  gr->Draw("AP");
  gr->SetMinimum(0);
  gr->Draw("AP");
  TLegend *myleg = mycan->BuildLegend(.7,.1,.9,.6);
  myleg->SetFillStyle(0);
  ofile->cd();
  gr->Write();
  ofile->cd();
  mycan->Write();
  if(draw != ""){
    mycan->SaveAs(string(draw + string(mycan->GetName())+".png").data());
  }
  delete mycan;
  
  get_eff_multigraphs(gr, scan_type, scan_nums, ofile, "simu", draw);
  if(scan_type == "Extraction"){
    get_eff_multigraphs(gr, scan_type, scan_nums, ofile, "Gushin", draw);
  }
  
  return;
}

void draw_hist_2d(TH2D* h, string scan_type, TFile *ofile, string directory){
  pair<TH2D*, TGraph*> norm_ratio = normalise_gain_histo_2d(h, scan_type);
  string field1 = scan_type.substr(0,scan_type.find_first_of("_"));
  string field2 = scan_type.substr(scan_type.find_first_of("_")+1);
  h->Draw("COLZ");
  h->GetZaxis()->SetRangeUser(0,3);
  h->SetTitle(string(string(h->GetName()) + ";" + field1 +";" + field2).data());
  gPad->Update();
  ofile->cd();
  h->Write();
  gPad->SaveAs(string(directory + string(h->GetName()) + ".png").data());
  gPad->SetName(string(string(h->GetName()) + "_pad").data());
  ofile->cd();
  gPad->Write();
  delete h;
  delete gPad;
  
  if(!norm_ratio.first){return;}
  norm_ratio.first->Draw("COLZ");
  norm_ratio.first->GetZaxis()->SetRangeUser(0,3);
  norm_ratio.first->SetTitle(string(string(norm_ratio.first->GetName()) + ";" + field1 +";" + field2).data());
  gPad->Update();
  ofile->cd();
  norm_ratio.first->Write();
  gPad->SaveAs(string(directory + string(norm_ratio.first->GetName()) + ".png").data());
  gPad->SetName(string(string(norm_ratio.first->GetName()) + "_pad").data());
  ofile->cd();
  gPad->Write();
  delete norm_ratio.first;
  delete gPad;
  
  norm_ratio.second->Draw("AP*");
  norm_ratio.second->SetMinimum(0);
  norm_ratio.second->GetXaxis()->SetRangeUser(0,norm_ratio.second->GetXaxis()->GetXmax());
  norm_ratio.second->Draw("AP*");
  gPad->Update();
  ofile->cd();
  norm_ratio.second->Write();
  gPad->SaveAs(string(directory + string(norm_ratio.second->GetName()) + ".png").data());
  gPad->SetName(string(string(norm_ratio.second->GetName()) + "_pad").data());
  ofile->cd();
  gPad->Write();
  delete norm_ratio.second;
  delete gPad;
  
  return;
}

void fit_gain(TGraphErrors* gr, vector<double> par){
  #if verbose
  cout << "  Fitting Gain vs LEM field for " << gr->GetName() << "..." << endl;
  #endif
  //fix rho, x and V, fit A, T and B. rho=PM/RT ->M and R constants, P and T from slow control (fix it at a constant value for all runs first, could apply correction later)
  //T ~ 1 Arho ~ 4000cm-1 Brho ~ 200kV/cm
  TF1* gain_fit = new TF1("gain_fit","[0]*exp([1]*0.1*exp(-[2]/x))",0,40);
  gain_fit->SetParameter(0,par[0]);
  gain_fit->SetParLimits(0,0,1);
  gain_fit->SetParameter(1,par[1]);
  gain_fit->SetParameter(2,par[2]);
  gr->Fit(gain_fit);
  delete gain_fit;
  return;
}

void fit_charging_up(TGraphErrors &gr){
  double begintime = 0;
  double endtime = 0;
  double begingain = 0;
  double endgain = 0;
  double a = 0;
  double b = 0;
  gr.GetPoint(0,begintime,begingain);
  gr.GetPoint(gr.GetN()-1,endtime,endgain);
  TF1* charging_up_fit = new TF1("charging_up_fit","[0]*x+[1]",0,endtime);
  a = (endgain-begingain)/(endtime-begintime);
  b = begingain-a*begintime;
  charging_up_fit->SetParameter(0,a);
  charging_up_fit->SetParameter(1,b);
//  TF1* charging_up_fit = new TF1("charging_up_fit","[0]/(1-(1-[0]/[1])*exp(-x*[2]))",0,endtime);
//  charging_up_fit->SetParameter(0,begingain*0.4);//from 3L paper (gain after charging up more or less 40% of initial gain)
//  charging_up_fit->SetParameter(1,begingain);
//  charging_up_fit->SetParameter(2,45000.);//from 3L paper (tau more or less half a day)
  gr.Fit(charging_up_fit);
  delete charging_up_fit;
  return;
}

double theoretical_gain(double T, double rho, double E){
  if(rho_ref == 0){
    cout << "gain_correction_for_rho: Please specify a reference density" << endl; return 0;
  }
  if(Arho_from_3L.first < 0){
    load_fit_3L();
  }
  double A = Arho_from_3L.first/rho_ref;
  double B = Brho_from_3L.first/rho_ref;
  return T*TMath::Exp(A*rho*0.1*TMath::Exp(-B*rho/E));
}

double gain_correction_for_rho(double rho, double E){
  double G0 = theoretical_gain(1,rho,E);
  double G0_corr = theoretical_gain(1,rho_ref,E);
  return G0_corr/G0;
}

double correct_dx_for_lifetime(double dx, double e_lifetime){
  return TMath::Exp(dx/(e_lifetime*drift_velocity));
}

double correct_for_drift(double E){
  return RFromBirk(E)/RFromBirk(0.5);
}

double get_density_for_hit(int hit_time, int run){
  if(Gr_Rho.GetN() == 0){
    if(!load_gr_and_h_rho(run)){return -1;}
  }
  return Gr_Rho.Eval(hit_time);
}

int read_or_do_fit(string filename_fitted, string filename_nonfitted, bool recreate_fit_file, string &ifilename, int &files_not_found){
  int i_read_fit = 0;
  if( ExistTest(filename_fitted) and !recreate_fit_file){
    struct stat t_stat_fitted;
    stat(filename_fitted.data(), &t_stat_fitted);
    struct stat t_stat_nonfitted;
    stat(filename_nonfitted.data(), &t_stat_nonfitted);
    if(t_stat_nonfitted.st_ctime < t_stat_fitted.st_ctime){
      #if verbose
      cout << "      Opening file: " << filename_fitted << endl;
      #endif
      ifilename = filename_fitted;
      i_read_fit = 1;
    }
    else{
      #if verbose
      cout << "      File fitted is older than file non fitted. Redoing fit." << endl;
      #endif
      ifilename = filename_nonfitted;
    }
  }
  else if( ExistTest(filename_nonfitted) ){
    #if verbose
    cout << "      Opening file: " << filename_nonfitted << endl;
    #endif
    ifilename = filename_nonfitted;
  }
  else{
    #if verbose
    cout << "    File " << filename_nonfitted << " not found." << endl;
    #endif
    files_not_found++;
    return 3;
  }
  return i_read_fit;
}


bool init_histo_track_cuts(vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected, int time){
  
  if( dx >= 100. ){
    #if verbose
    cout << "ERROR: can not have dx superior to total height of detector" << endl;
    #endif
    return false;
  }
  string title = "dQds";
  if(time > 0){title = to_string(time);}
  string histname = "dQds_view0";
  hdQds.push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
  histname = "dQds_view1";
  hdQds.push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
  histname = "dQds_view0_Dx_Corrected";
  hdQds_Dx_Corrected.push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
  histname = "dQds_view1_Dx_Corrected";
  hdQds_Dx_Corrected.push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
  
  for( int X = 0; X < (int)(100/dx); X++ ){
    int x = dx*X;
    for( auto lem : lems ){
      if(x == 0){
        histname = "dQds_LEM_"+to_string(lem)+"_view0";
        hdQds_ByLEMs[lem].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
        histname = "dQds_LEM_"+to_string(lem)+"_view1";
        hdQds_ByLEMs[lem].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
        histname = "dQds_LEM_"+to_string(lem)+"_view0_Dx_Corrected";
        hdQds_ByLEMs_Dx_Corrected[lem].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
        histname = "dQds_LEM_"+to_string(lem)+"_view1_Dx_Corrected";
        hdQds_ByLEMs_Dx_Corrected[lem].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
      }
      histname = "dQds_dx_"+to_string(x)+"_LEM_"+to_string(lem)+"_view0";
      hdQds_ByDx_ByLEMs[x][lem].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
      histname = "dQds_dx_"+to_string(x)+"_LEM_"+to_string(lem)+"_view1";
      hdQds_ByDx_ByLEMs[x][lem].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
    }
    histname = "dQds_dx_"+to_string(x)+"_view0";
    hdQds_ByDx[x].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
    histname = "dQds_dx_"+to_string(x)+"_view1";
    hdQds_ByDx[x].push_back(TH1D(histname.data(), title.data(), 100, 0, 50));
  }
  return true;
}

void rec_track_dQds(track t, vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected){
  for(auto h : t.hits_trk){
    if( !isGood_lem(h.lem) ){continue;}
    if( !IsGoodChan(h) ){continue;}
    double dq;
    double ds;
    if(method_ds == "3D"){ds = h.dx_3D;}
    else if(method_ds == "local"){ds = h.dx_local;}
    else{ds = h.dx_phil;}
    if(method_dQ == "sum"){dq = h.dq_sum;}
    else{dq = h.dq_integral;}
//    double dqds = h.gain_density_correction_factor * dq/ds;
    double dqds = dq/ds;
    
    //cut on dqdx
    if( h.sp_x < tpc_boundaries[0] or h.sp_x >= tpc_boundaries[1] or h.sp_y < tpc_boundaries[2] or h.sp_y > tpc_boundaries[3] or h.sp_z < tpc_boundaries[4] or h.sp_z > tpc_boundaries[5] ){ continue; }
    
    
    hdQds[h.view].Fill(dqds);
    hdQds_ByLEMs[h.lem][h.view].Fill(dqds);
    hdQds_Dx_Corrected[h.view].Fill(dqds*h.purity_correction);
    hdQds_ByLEMs_Dx_Corrected[h.lem][h.view].Fill(dqds*h.purity_correction);
    hdQds_ByDx[(int)((int)((h.sp_x+50)/dx)*dx)][h.view].Fill(dqds);
    hdQds_ByDx_ByLEMs[(int)((int)((h.sp_x+50)/dx)*dx)][h.lem][h.view].Fill(dqds);
  }//end hits;
  return;
}

bool select_tracks(string cut_type, vector<track> tracks, vector<track> & mips, vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected,  map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected){
  //select particles crossing the detector in any direction
  int count_mip=0;
  bool reverse_cut = false;
  
  if(cut_type.find("reverse") != string::npos){
    reverse_cut = true;
  }
  
  for(auto t : tracks){

    double mag = t.length;
    
    if(cut_type == "before_cuts"){
      rec_track_dQds(t, hdQds, hdQds_ByLEMs, hdQds_ByDx, hdQds_ByDx_ByLEMs, hdQds_Dx_Corrected, hdQds_ByLEMs_Dx_Corrected);
      mips.push_back(t);
      count_mip++;
      continue;
    }
    
    if( (!reverse_cut and
        mag > length_cut
        and t.id < ntracks
        and t.event_ntracks < nmaxtracksinevent
        and t.nhitsview[0] > nhits_cut
        and t.nhitsview[1] > nhits_cut
        and (t.throughgoing or !only_throughgoing)
        and (t.on_edge or !keep_on_edge)
        and (t.throughgoingY or !only_throughgoingY)
        and (t.on_edgeX or !keep_on_edgeY)
        and (t.throughgoingZ or !only_throughgoingZ)
        and (t.on_edgeX or !keep_on_edgeZ)
        and (t.throughgoingX or !only_throughgoingX)
        and (t.on_edgeX or !keep_on_edgeX)
        and (t.highway_selected or !highway)
        and (t.blc or !below_central_lem)
        and (t.theta < -theta_cut or theta_cut >= 0)
        and (t.theta > theta_cut or theta_cut <= 0)
        and get_reduced_phi(t) > phi_range[0]
        and get_reduced_phi(t) < phi_range[1]
        and t.theta > theta_range[0]
        and t.theta < theta_range[1]
        and (fabs(t.start_theta - t.end_theta) < deltatheta_cut)
        and (t.long_in_view >= 0 or !only_long_view)
        and (fabs(t.start_phi - t.end_phi) < deltaphi_cut) )
        
        or (reverse_cut and
        !(  mag > length_cut
        and t.id < ntracks
        and t.event_ntracks < nmaxtracksinevent
        and t.nhitsview[0] > nhits_cut
        and t.nhitsview[1] > nhits_cut
        and (t.throughgoing or !only_throughgoing)
        and (t.on_edge or !keep_on_edge)
        and (t.throughgoingY or !only_throughgoingY)
        and (t.on_edgeX or !keep_on_edgeY)
        and (t.throughgoingZ or !only_throughgoingZ)
        and (t.on_edgeX or !keep_on_edgeZ)
        and (t.throughgoingX or !only_throughgoingX)
        and (t.on_edgeX or !keep_on_edgeX)
        and (t.highway_selected or !highway)
        and (t.blc or !below_central_lem) )
        and (t.theta < -theta_cut or theta_cut >= 0)
        and (t.theta > theta_cut or theta_cut <= 0)
        and get_reduced_phi(t) > phi_range[0]
        and get_reduced_phi(t) < phi_range[1]
        and t.theta > theta_range[0]
        and t.theta < theta_range[1]
        and (fabs(t.start_theta - t.end_theta) < deltatheta_cut)
        and (t.long_in_view >= 0 or !only_long_view)
        and (fabs(t.start_phi - t.end_phi) < deltaphi_cut) )
        ) {
      
      //cut on the track angle. Avoid track parallel to a view, or parallel to drift direction. Also ignore bended tracks 
//      if( theta_cut > 0 and (t.theta < theta_cut or t.theta > 180-theta_cut) ){continue;}
//      if( phi_cut > 0 and (fabs(t.phi) - ((int)(fabs(t.phi)/90))*90 < phi_cut or fabs(t.phi) - ((int)(fabs(t.phi)/90))*90 > 90-phi_cut) ){continue;}
      bool Continue = false;
      t.nhitsview[0] = 0; t.nhitsview[1] = 0;
      t.sumQ[0] = 0; t.sumQ[1] = 0;
      
      //remove hits with too small ds
      for( vector<hit>::iterator h = t.hits_trk.begin(); h != t.hits_trk.end(); ){
        double ds;
        double dQ;
        if(method_dQ == "sum"){dQ = h->dq_sum;}
        else{dQ = h->dq_integral;}
        if(method_ds == "3D"){ds = h->dx_3D;}
        else if(method_ds == "local"){ds = h->dx_local;}
        else{ds = h->dx_phil;}
        double dqds = dQ/ds;
        if(
        (!reverse_cut and ( dqds <= dQdx_cut_min or dqds > dQdx_cut_max or ds > ds_cut or ds < pitch or h->sp_x < MinX or h->sp_x > MaxX or h->sp_y < MinY or h->sp_y > MaxY or h->sp_z < MinZ or h->sp_z > MaxZ or h->GoF > GoF_cut or !isGood_lem(h->lem) or !IsGoodChan(h->sp_y, h->sp_z) or (only_long_view and t.long_in_view != h->view) ) ) 
        or 
        (reverse_cut and !(dqds <= dQdx_cut_min or dqds > dQdx_cut_max or ds > ds_cut or ds < pitch or h->sp_x < MinX or h->sp_x > MaxX or h->sp_y < MinY or h->sp_y > MaxY or h->sp_z < MinZ or h->sp_z > MaxZ or h->GoF > GoF_cut or !isGood_lem(h->lem) or !IsGoodChan(h->sp_y, h->sp_z) or (only_long_view and t.long_in_view != h->view) ) ) 
        ){
          h = t.hits_trk.erase(h);
        }
        else{
          t.nhitsview[h->view]++;
          t.sumQ[h->view] += dQ;
          ++h;
        }
      }
      if(t.hits_trk.size() == 0){continue;}
      double Qasym = 0;
      if(t.sumQ[0] > 0){
        Qasym = t.sumQ[1]/t.sumQ[0];
      }
      if( (!reverse_cut and (Qasym < Qasym_cut[0]       or Qasym > Qasym_cut[1]      )) or (reverse_cut and !(Qasym < Qasym_cut[0]       or Qasym > Qasym_cut[1]      )) ){continue;}
      if( (!reverse_cut and (t.nhitsview[0] < nhits_cut or t.nhitsview[1] < nhits_cut)) or (reverse_cut and !(t.nhitsview[0] < nhits_cut or t.nhitsview[1] < nhits_cut)) ){continue;}
      mips.push_back(t);
      rec_track_dQds(t, hdQds, hdQds_ByLEMs, hdQds_ByDx, hdQds_ByDx_ByLEMs, hdQds_Dx_Corrected, hdQds_ByLEMs_Dx_Corrected);
      count_mip++;
    }//end mag
  }//end for tracks
  
  #if verbose
  cout << " Selected " << count_mip << " tracks over " << tracks.size() << " tracks " <<  endl;
  #endif
  return true;
}

//check if all the directoryies in a file address exist, creates them if possible. Can take either a file with complete address or a directory
bool check_and_mkdir(string path){
  istringstream dummy(path);
  vector<string> subpaths((istream_iterator<WordDelimitedBySlash>(dummy)),istream_iterator<WordDelimitedBySlash>());
  string partpath = "";
  unsigned int i = 0;
  for(auto sp : subpaths){
    string previous = partpath;
    ++i;
    if(partpath == "" and sp == ""){
      partpath += "/";
      continue;
    }
    else if(partpath == ""){return true;}
    else{partpath += sp;}
    if(partpath.find(".") != string::npos and i == subpaths.size()){return true;}
  
    if(!ExistTest(partpath)){
      if(access(previous.data(),W_OK) != 0){
        cout << "ERROR in check_and_mkdir: do not have permission to create " << partpath << " directory." << endl;
        return false;
      }
      #if verbose
      cout << "Creating directory " << partpath << endl;
      #endif
      mkdir(partpath.data(),0777);
    }
    partpath += "/";
  }
  return true;
}

bool get_histo_in_inputfile(TH1D &hdQds, TFile *runfile, string name_to_get, bool &read_fit){
  TH1D *dummyh = 0;
  runfile->GetObject(name_to_get.data(), dummyh);
  if(dummyh == 0){
    #if verbose
    cout << "    Could not find " << name_to_get << " in " << runfile->GetName() << endl;
    #endif
    if(string(runfile->GetName()).find("fitted") != string::npos){
      #if verbose
      cout << "    Will try using non-fitted file instead" << endl;
      #endif
      read_fit = false;
      return true;
    }
    else{
      cout << "Please check your dQds file" << endl;
      return false;
    }
  }
  hdQds = TH1D(*dummyh);
  hdQds.GetListOfFunctions()->Clear();
  if(dummyh->GetListOfFunctions()->GetSize() > 0){
    hdQds.GetListOfFunctions()->Add(dummyh->GetListOfFunctions()->At(0));
  }
  return true;
}

string set_cuts(string cut_type){
  string to_return = "";
  if(cut_type.find("reverse") != string::npos){
    to_return = "reverse_";
  }
  if(cut_type.find("laura") != string::npos){
    length_cut = 50.;
    nmaxtracksinevent = 80;
    only_throughgoingX = true;
    keep_on_edgeX = true;
    to_return = to_return + "laura_";
    
  }
  else if(cut_type.find("christoph") != string::npos){
    nmaxtracksinevent = 10;
    only_throughgoingX = true;
    below_central_lem = true;
    to_return = to_return + "christoph_";
    
  }
  else if(cut_type.find("philippe") != string::npos){
    nmaxtracksinevent = 10;
    only_throughgoingX = true;
//    keep_on_edgeX = true;
    theta_cut = 2.3;
    to_return = to_return + "philippe_";
    
  }
  if(cut_type.find("dray") != string::npos){
    dray_miti = true; to_return+="dray_";
  }
  if(cut_type.find("common") != string::npos){
    length_cut = 50.;method_dQ = "sum";method_ds = "local"; to_return+="common_"; highway = true;
  }
  if(cut_type.find("length") != string::npos and cut_type.find("laura") == string::npos) {
    length_cut = 50.; to_return+="length_";
  }
  if(cut_type.find("nolen") != string::npos) {
    length_cut = 0.; to_return+="nolen_";
  }
  if(cut_type.find("Ds") != string::npos) {
    ds_cut = 1.; to_return+="Ds_";
  }
//  if(cut_type.find("GoF") != string::npos) {GoF_cut = 0.2; to_return+="GoF_";}
  if(cut_type.find("tgx") != string::npos) {
    only_throughgoingX = true; to_return+="tgx_";
  }
  if(cut_type.find("tgy") != string::npos) {
    only_throughgoingY = true; to_return+="tgy_";
  }
  if(cut_type.find("tgz") != string::npos) {
    only_throughgoingZ = true; to_return+="tgz_";
  }
  if(cut_type.find("tg") != string::npos and cut_type.find("tgx") == string::npos and cut_type.find("tgy") == string::npos and cut_type.find("tgz") == string::npos) {
    only_throughgoing = true; to_return+="tg_";
  }
  if(cut_type.find("on_edgeX") != string::npos) {
    keep_on_edgeX = true; to_return+="on_edgeX_";
  }
  if(cut_type.find("on_edgeY") != string::npos) {
    keep_on_edgeY = true; to_return+="on_edgeY_";
  }
  if(cut_type.find("on_edgeZ") != string::npos) {
    keep_on_edgeZ = true; to_return+="on_edgeZ_";
  }
  if(cut_type.find("on_edge") != string::npos and cut_type.find("on_edgeX") == string::npos and cut_type.find("on_edgeY") == string::npos and cut_type.find("on_edgeZ") == string::npos) {
    keep_on_edge = true; to_return+="on_edge_";
  }
  if(cut_type.find("below_central_lem") != string::npos or cut_type.find("bcl") != string::npos) {
    below_central_lem = true; to_return+="bcl_";
  }
  if(cut_type.find("ntracks") != string::npos) {
    ntracks = 10; to_return+="ntracks_";
  }
  if(cut_type.find("maxtracksevent") != string::npos) {
    nmaxtracksinevent = 10; to_return+="maxtracksevent_";
  }
  if(cut_type.find("mtheta") != string::npos) {
    theta_cut = -2.3; to_return+="mtheta_";
  }
  if(cut_type.find("ptheta") != string::npos) {
    theta_cut = 2.3; to_return+="ptheta_";
  }
  if(cut_type.find("deltatheta") != string::npos) {
    deltatheta_cut = 10.; to_return+="deltatheta_";
  }
  if(cut_type.find("deltaphi") != string::npos) {
    deltaphi_cut = 10.; to_return+="deltaphi_";
  }
  if(cut_type.find("phirange") != string::npos) {
    phi_range[0] = 30; phi_range[1] = 60; to_return+="phirange_";
  }
  if(cut_type.find("thetarange") != string::npos) {
    theta_range[0] = 135; theta_range[1] = 180; to_return+="thetarange_";
  }
  if(cut_type.find("Qasym") != string::npos) {
    Qasym_cut[0] = 0.8; Qasym_cut[1] = 1.2; to_return+="Qasym_";
  }
  if(cut_type.find("nhits") != string::npos) {
    nhits_cut = 50; to_return+="nhits_";
  }
  if(cut_type.find("onlylongview") != string::npos) {
    only_long_view = true; to_return+="onlylongview_";
  }
  
  if(to_return == ""){to_return = "before_cuts";}
  else{to_return = to_return.substr(0, to_return.size()-1);}
  return to_return;
}

void set_bad_runs(){
  bad_runs["Ds_local_sum"] = {771,992,1005,1165,1166,1167,1172,1173,1184,1185,1195};
  bad_runs_lems["Ds_local_sum"][1002] = {2};
  bad_runs_lems["Ds_local_sum"][1037] = {6};
  bad_runs_lems["Ds_local_sum"][1039] = {2,4};
  bad_runs_lems["Ds_local_sum"][1174] = {9,11};
  bad_runs_lems["Ds_local_sum"][1175] = {9,11};
  bad_runs_lems["Ds_local_sum"][1177] = {9,11};
  bad_runs_lems["Ds_local_sum"][1183] = {2,9};
  bad_runs_lems["Ds_local_sum"][1187] = {9,11};
  bad_runs_lems["Ds_local_sum"][1188] = {2,9,11};

  bad_runs_lems["common_local_sum"][783] = {7,8};
  bad_runs_lems["common_local_sum"][786] = {8};
  bad_runs_lems["common_local_sum"][785] = {4};
  bad_runs_lems["common_local_sum"][784] = {4,5,8};
  bad_runs_lems["common_local_sum"][840] = {5,7,8};
  bad_runs_lems["common_local_sum"][835] = {8};
  bad_runs_lems["common_local_sum"][833] = {8};
  
  bad_runs_yz["zob"][1] = {make_pair(-1,-1)};

  return;
}

bool pressure(string srun){
  
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
  if(!ExistTest(path)){cout << "ERROR: file " << path << " not found" << endl; return false;}
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
    return false;
  }
  int i = 0;
  int npar = 0;
  vector<int> columns_to_rec = {};
  map<string,vector<double> > map_par_values;
  vector<double> mean_cryo_temps = {};
  vector<double> cryo_pos_for_plot = {};
  double rho_var = 0;
  double rho_mean = 0;
  
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
      npar = col;
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
    if(npar != col){
      cout << "ERROR in pressure: run has incoherent number of parameters: " << col << " instead of " << npar << endl;
      return false;
    }
  }
  ifile.close();
  
  if(map_par_values["date"].size() == 0){
    cout << "Error in pressure: empty run" << endl;
    return false;
  }
  
  cout << "...done" << endl;
  TFile ofile_pressure(string(directory + "pressure.root").data(),"RECREATE");
  TFile ofile_temperature_CRP(string(directory + "temperature_CRP.root").data(),"RECREATE");
  TFile ofile_temperature_cryostat(string(directory + "temperature_cryostat.root").data(),"RECREATE");
  
  
  for(auto par : map_par_values){
////////////////////////////////////////////////
    //param vs time
////////////////////////////////////////////////

    cout << "Plotting " << par.first << "..." << endl;
    if(par.first == "date"){continue;}
    TGraph gr(par.second.size(), &map_par_values["date"][0], &par.second[0]);
    gr.SetTitle(string(par.first).data());
    gr.SetName(string("graph_"+par.first).data());
    gr.Draw("AC");
    gr.GetXaxis()->SetTimeDisplay(1);
    gr.GetXaxis()->LabelsOption("v");
    if(map_par_values["date"].back()-map_par_values["date"][0] > 24*3600){
      gr.GetXaxis()->SetTimeFormat("%d/%m %Hh%M");
      gr.GetXaxis()->SetNdivisions(5);
    }
    else{
      gr.GetXaxis()->SetTimeFormat("%Hh%M");
      gr.GetXaxis()->SetNdivisions(8);
    }
    gr.Draw("AC");
    gPad->Update();
    gPad->Modified();
    if(par.first == "PE0006"){ofile_pressure.cd();}
    else if(find(temp_crp.begin(), temp_crp.end(), par.first) != temp_crp.end()){ofile_temperature_CRP.cd();}
    else if(temp_cryostat.find(par.first) != temp_cryostat.end()){ofile_temperature_cryostat.cd();}
    gr.Write();
    gPad->SetName(string(string(gr.GetName()) + "_pad").data());
    gPad->Write();
//    gPad->SaveAs(string(directory + "graph_"+par.first + ".png").data());
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
    gPad->Write();
//    gPad->SaveAs(string(directory + "histo_"+par.first + ".png").data());
    delete gPad;
    
    if(par.first != "PE0006"){
    
////////////////////////////////////////////////
    //density vs time
////////////////////////////////////////////////

      vector<double> rho = {};
      for(unsigned int i = 0; i < par.second.size(); i++){
        rho.push_back(map_par_values["PE0006"][i]/par.second[i]);
      }
      TGraph gr_rho(rho.size(), &map_par_values["date"][0], &rho[0]);
      gr_rho.SetTitle(string("PE0006/"+par.first).data());
      gr_rho.SetName(string("graph_PE0006Over"+par.first).data());
      gr_rho.Draw("AC");
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
      gPad->Write();
//      gPad->SaveAs(string(directory + "graph_PE0006Over"+par.first + ".png").data());
      delete gPad;
      
////////////////////////////////////////////////
    //density distribution
////////////////////////////////////////////////
      
      TH1D h_rho(string("histo_PE0006Over"+par.first).data(),string("PE0006Over"+par.first+";PE0006Over"+par.first+";#").data(),100,0.99*rho[min_element(rho.begin(),rho.end()) - rho.begin()],rho[max_element(rho.begin(),rho.end()) - rho.begin()]*1.01);
      for(auto value : rho){
        h_rho.Fill(value);
      }
      if(par.first == "TE0056"){
        rho_mean = h.GetMean();
        rho_var = h.GetStdDev();
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
    gains_middle.push_back(theoretical_gain(1,rho_mean,fiel));
    gains_up.push_back(theoretical_gain(1,rho_mean + rho_var/2.,fiel));
    gains_down.push_back(theoretical_gain(1,rho_mean - rho_var/2.,fiel));
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
  gain_up.SetTitle(string("gain with density variating from " + to_string((int)rho_mean - rho_var/2.) + " to " + to_string((int)rho_mean + rho_var/2.)+";LEM Field (kV/cm);Gain in LEM").data());
  gain_up.Draw("AC");
  gain_middle.Draw("SAME C");
  gain_down.Draw("SAME C");
  gPad->SetName("gain_variations");
  gPad->SaveAs(string(directory + "gain_variations.png").data());
  ofile_temperature_cryostat.cd();
  gPad->Write();
  delete gPad;
      
  ofile_temperature_cryostat.Close();
  
  return true ;
}

bool save_run_header(TMyFileHeader header){
  string ofilename = runs_headers+to_string(header.GetRun())+".root";
  TFile ofile(ofilename.data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
    return false;
  }
  #if verbose
  cout << " Header header_" << to_string(header.GetRun()).data() << " saved in " << ofilename << "." << endl;
  #endif
  
  header.Write(string("header_"+to_string(header.GetRun())).data());
  ofile.Close();
  return true;
}

bool set_gain_processed(vector<int> run_list, string cut_type, string m_dq, string m_ds){
  method_dQ = m_dq;
  method_ds = m_ds;
  set_bad_runs();
  if(!Load_Version("July")){return false;}
  if (run_list.size() == 0){
    #if verbose
    cout << "Processing all runs in " << dQds_Output << "*..." << endl;
    #endif
    string wildcard_path = dQds_Output + cut_type + "/*";
    for( auto irun : glob(wildcard_path) ){
      if(irun.find("plots") != string::npos){continue;}
      if(irun.find("fitted") != string::npos){continue;}
      irun = irun.substr(0, irun.find_first_of("."));
      run_list.push_back(atoi(irun.data()));
    }
  }
  string cut_type_and_methods = cut_type + "_" + method_ds + "_" + method_dQ;
  for(auto run : run_list){
    bool good = IsGood_run(cut_type_and_methods, run);
    string ofilename = runs_headers+to_string(run)+".root";
    TFile ofile(ofilename.data(), "UPDATE");
    if(!ofile.IsOpen()){
      cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
      return false;
    }
    string headername = "header_"+to_string(run);
    TMyFileHeader *h = 0;
    ofile.GetObject(headername.data(),h);
    if(h->AreGainsProcessed().find(cut_type_and_methods) == h->AreGainsProcessed().end()){cout <<"No cut named " << cut_type_and_methods << " found." << endl; return false;}
    h->SetGainProcessed(cut_type_and_methods, good);
    for(auto lem : lems){
      bool goodlem = IsGood_lem_gain(cut_type_and_methods, run, lem);
      h->SetGainProcessedLEM(cut_type_and_methods, lem, good and goodlem);
    }
    if(h->GetYZMPVs().size() != 0){
      vector<int> begins = {-(int)(lem_size/dy), 0};
      vector<int> ends = {(int)(lem_size/dy), (int)(6*lem_size/dz)};
      for( int Y = begins[0]; Y < ends[0]; Y++ ){
        for( int Z = begins[1]; Z < ends[1]; Z++ ){
          pair<int,int> YZ = make_pair(Y,Z);
          double mpv0 = h->GetYZMPVs()[cut_type_and_methods][YZ].first;
          double mpv1 = h->GetYZMPVs()[cut_type_and_methods][YZ].second;
          bool goodYZ = (IsGood_yz(cut_type_and_methods, run, make_pair(Y,Z)) and mpv0 > 0 and mpv1 > 0);
          h->SetYZGainProcessed(cut_type_and_methods, make_pair(Y,Z), goodYZ);
        }
      }
    }
    h->Write(headername.data(),TObject::kOverwrite);
    #if verbose
    cout << " Header " << headername << " updated in " << ofilename << "." << endl;
    #endif
    ofile.Close();
  }
  return true;
}

bool save_runs_headers(vector<int> run_list){
  if(version == ""){if(!Load_Version("Feb")){return false;}}
  if(runs_and_fields.size() == 0){if(!load_run_lists()){return false;}}
  
  string path = path_wa105_311data;
  string path_output = runs_headers;
  if (run_list.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path_wa105_311data << "*..." << endl;
    #endif
    string wildcard_path = path + "*";
    for( auto irun : glob(wildcard_path) ){
      run_list.push_back(atoi(irun.data()));
    }
  }
  string filename = "";
  unsigned int failed = 0;
  
  for(auto run : run_list){
  
    string ofilename = path_output+to_string(run)+".root";
    double drift = runs_and_fields[run]["Drift"];
    double amplification = runs_and_fields[run]["Amplification"];
    double extraction = runs_and_fields[run]["Extraction"];
    double induction = runs_and_fields[run]["Induction"];
    double extr_eff_simu = get_eff(extraction, amplification);
    double extr_eff_gushin = get_eff(extraction);
    double ind_eff = get_eff(-1,amplification,induction);
    double drift_correction_factor = correct_for_drift(drift);
    string trigger = runs_triggers[run];
    int tstart = 0;
    int tend = 0;
    
    TChain chain("analysistree/anatree");
    string foldername = path+to_string(run)+"/";
    int n_subruns = glob(string(foldername+"*").data()).size();
    if( n_subruns == 0){
      #if verbose
      cout << foldername << " not found or empty" << endl;
      #endif
      failed++;
      continue;
    }
    #if verbose
    cout << "Processing " << n_subruns << " subruns in run " << run << endl;
    #endif

    string file = "";
    for(int i=0; i< n_subruns; i++){
      
      file=to_string(run)+"-"+to_string(i);
      if(file == "1009-21"){failed++;continue;}
      if(path_311data.find("Feb") != string::npos){
        file = file + "-Parser.root";
      }
      else if(path_311data.find("June") != string::npos){
        file = file + "-RecoFast-Parser.root";
      }
      else if(path_311data.find("July") != string::npos){
        file = file + "-RecoFast-Parser.root";
      }
      else{
        cout << "ERROR: unknown reconstruction version" << endl;
        return false;
      }
      filename=foldername+file;
      if( ExistTest(filename) ){
        #if verbose
        cout << "Adding file: " << filename << endl;
        #endif
        chain.Add(filename.data());
      }
      else{
        #if verbose
        cout << "File " << filename << " not found." << endl;
        #endif
        failed++;
        continue;
      }
    }
    #if verbose
    cout << "Files added. " << endl;
    #endif
    
    int NEntries = (int)chain.GetEntries();
    int  tEventTimeSeconds;
    chain.SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
    chain.GetEntry(0);
    tstart = tEventTimeSeconds;
    chain.GetEntry(NEntries-1);
    tend = tEventTimeSeconds;
    if(!ExistTest(slow_control + to_string(run) + "/" + to_string(run) + ".txt") and !IsBatch){
      #if verbose
      cout << "Getting slow control data for run  " << run << "..." << endl;
      #endif
      string command  = path_311analysis+"slow_control/Pcotte_getdata.sh " + to_string(tstart) + " " + to_string(tend) + " " + to_string(run);
      system(command.data());
    }
    if(!ExistTest(slow_control + "/" + to_string(run) + "/" + to_string(run) + ".txt")){
      cout << "ERROR: can not get slow control data for run " << run << endl;
      failed++;
      continue;
    }
    #if verbose
    cout << "Opening slow control data for run  " << run << "..." << endl;
    #endif
    if(Gr_Rho.GetN() == 0 or H_Rho.GetEntries() == 0){
      if(!load_gr_and_h_rho(run)){failed++;continue;}
    }
    double rho = H_Rho.GetMean();
    double rho_var = H_Rho.GetStdDev();
    double density_correction_factor = gain_correction_for_rho(rho, amplification);
    
    TMyFileHeader header(run, drift, extraction, amplification, induction, tstart ,tend,string("header_"+to_string(run)).data(), rho, rho_var, extr_eff_simu, 0, extr_eff_gushin, 0, ind_eff, 0, drift_correction_factor, 0, density_correction_factor, 0, theoretical_gain(1,rho,amplification), theoretical_gain(1,rho_ref,amplification), trigger);
    TFile ofile(ofilename.data(), "RECREATE");
    if(!ofile.IsOpen()){
      cout << " ERROR: TFile " << ofilename << " can't be created " << endl;
      failed++;
      continue;
    }
    #if verbose
    cout << " Header header_" << to_string(run).data() << " saved in " << ofilename << " has been created " << endl;
    #endif
    
    header.Write(string("header_"+to_string(run)).data());
    ofile.Close();
    
  }//for run
  if(failed == run_list.size()){
    cout << "ERROR in save_run_headers: all runs failed" << endl;
    return false;
  }
  return true;
}


