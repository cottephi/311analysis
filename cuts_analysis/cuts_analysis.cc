#include <311Lib.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

void analyse_cuts(vector<track> tracks, string outpath){

  int run = tracks[0].run;
  string dirname = outpath;
  outpath = outpath + "analysis.root";
  TFile ofile(outpath.data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << "ERROR: can not open " << outpath << endl;
    return;
  }
  
  map<int,int> events_and_number_of_tracks;
  
  int nb_tracks = 0;
  int nb_hits_view0 = 0;
  int nb_hits_view1 = 0;
  
  int max_hits_per_event = 0;
  int max_tracks_per_event = 0;
  
  TH1D* number_of_hits_per_tracks_view0 = new TH1D("number_of_hits_per_tracks_view0","tracks_and_number_of_hits_per_tracks_view0",100,0,1120);
  TH1D* number_of_hits_per_tracks_view1 = new TH1D("number_of_hits_per_tracks_view1","tracks_and_number_of_hits_per_tracks_view1",100,0,1120);
  TH1D* dQds_run_view0 = new TH1D("dQds_run_view0","dQds_run_view0",100,0,50);
  TH1D* dQds_run_view1 = new TH1D("dQds_run_view1","dQds_run_view1",100,0,50);
  TH1D* length_of_tracks = new TH1D("length_of_tracks","length_of_tracks",100,0,332);
  TH1D* theta_of_tracks = new TH1D("theta_of_tracks","theta_of_tracks",100,0,TMath::Pi());
  TH1D* phi_of_tracks = new TH1D("phi_of_tracks","phi_of_tracks",100,-TMath::Pi(),TMath::Pi());
  TH1D* ds_of_hits_view0 = new TH1D("ds_of_hits_view0","ds_of_hits_view0",100,0,6);
  TH1D* ds_of_hits_view1 = new TH1D("ds_of_hits_view1","ds_of_hits_view1",100,0,6);
  TH1D* GoF_of_hits_view0 = new TH1D("GoF_of_hits_view0","GoF_of_hits_view0",100,0,10);
  TH1D* GoF_of_hits_view1 = new TH1D("GoF_of_hits_view1","GoF_of_hits_view1",100,0,10);
  TH2D* xz = new TH2D("xz","xz",300, 0, 300,100,-50,50);
  TH2D* xy = new TH2D("xy","xy",100,-50,50, 100,-50,50);
  TH2D* yz = new TH2D("yz","yz",300, 0, 300,100,-50,50);
  TH3D* xyz = new TH3D("xyz","xyz",300,0,300,100,-50,50,100,-50,50);
  TH1D* ds_of_tracks_view0 = new TH1D("ds_of_tracks_view0","ds_of_tracks_view0",100,0,10);
  TH1D* ds_of_tracks_view1 = new TH1D("ds_of_tracks_view1","ds_of_tracks_view1",100,0,10);
  
  double ds;
  double dq;
  double dqds;
  
  for( auto t : tracks){
    nb_tracks++;
    int hits_track_view0 = 0;
    int hits_track_view1 = 0;
    if(events_and_number_of_tracks.find(t.event) == events_and_number_of_tracks.end()){
      events_and_number_of_tracks[t.event] = 1;
    }
    else{
      events_and_number_of_tracks[t.event]++;
      if( events_and_number_of_tracks[t.event] > max_tracks_per_event){
        max_tracks_per_event = events_and_number_of_tracks[t.event];
      }
    }
    for(auto h : t.hits_trk){

      if(method_ds == "3D"){ds = h.dx_3D;}
      else{ds = h.dx_local;}
      if(method_dQ == "sum"){dq = h.dq_sum;}
      else{dq = h.dq_integral;}
      dqds = h.gain_density_correction_factor * dq/ds;
      if(!isGood_lem(h.lem)){
        continue;
      }
      if(!IsGoodChan(h)){
        continue;
      }
      if(h.view == 0){
        dQds_run_view0->Fill(dqds);
        ds_of_hits_view0->Fill(ds);
        GoF_of_hits_view0->Fill(h.GoF);
        hits_track_view0++;
        nb_hits_view0++;
      }
      if(h.view == 1){
        dQds_run_view1->Fill(dqds);
        ds_of_hits_view1->Fill(ds);
        GoF_of_hits_view1->Fill(h.GoF);
        hits_track_view1++;
        nb_hits_view1++;
      }
    }
    xz->Fill(t.start_z,t.start_x);
    xz->Fill(t.end_z,t.end_x);
    xy->Fill(t.start_y,t.start_x);
    xy->Fill(t.end_y,t.end_x);
    yz->Fill(t.start_z,t.start_y);
    yz->Fill(t.end_z,t.end_y);
    xyz->Fill(t.start_z,t.start_y,t.start_x);
    xyz->Fill(t.end_z,t.end_y,t.end_x);
    number_of_hits_per_tracks_view0->Fill(hits_track_view0);
    number_of_hits_per_tracks_view1->Fill(hits_track_view1);
    length_of_tracks->Fill(t.length);
    theta_of_tracks->Fill(t.theta);
    phi_of_tracks->Fill(t.phi);
    ds_of_tracks_view0->Fill(get_dss(t)[0]);
    ds_of_tracks_view1->Fill(get_dss(t)[1]);
  }
  
  ///////////////////////////////////////save///////////////////////////////////////////
  
  TH1D* number_of_tracks_per_event = new TH1D("number_of_tracks_per_event","number_of_tracks_per_event",100,0,max_tracks_per_event);
  
  for(auto tracks : events_and_number_of_tracks){
    number_of_tracks_per_event->Fill(tracks.second);
  }
  ds_of_tracks_view0->Write();
  ds_of_tracks_view0->Draw();
  gPad->SaveAs(string(dirname + string(ds_of_tracks_view0->GetName())+".png").data());
  ds_of_tracks_view1->Write();
  ds_of_tracks_view1->Draw();
  gPad->SaveAs(string(dirname + string(ds_of_tracks_view1->GetName())+".png").data());
  xyz->Write();
  xyz->Draw("COL");
  gPad->SaveAs(string(dirname + string(xyz->GetName())+".png").data());
  xz->Write();
  xz->Draw("COL");
  gPad->SaveAs(string(dirname + string(xz->GetName())+".png").data());
  theta_of_tracks->Write();
  xy->Write();
  xy->Draw("COL");
  gPad->SaveAs(string(dirname + string(xy->GetName())+".png").data());
  yz->Write();
  yz->Draw("COL");
  gPad->SaveAs(string(dirname + string(yz->GetName())+".png").data());
  ds_of_hits_view0->Write();
  ds_of_hits_view0->Draw();
  gPad->SaveAs(string(dirname + string(ds_of_hits_view0->GetName())+".png").data());
  ds_of_hits_view1->Write();
  ds_of_hits_view1->Draw();
  gPad->SaveAs(string(dirname + string(ds_of_hits_view1->GetName())+".png").data());
  GoF_of_hits_view0->Write();
  GoF_of_hits_view0->Draw();
  gPad->SaveAs(string(dirname + string(GoF_of_hits_view0->GetName())+".png").data());
  GoF_of_hits_view1->Write();
  GoF_of_hits_view1->Draw();
  gPad->SaveAs(string(dirname + string(GoF_of_hits_view1->GetName())+".png").data());
  theta_of_tracks->Write();
  theta_of_tracks->Draw();
  gPad->SaveAs(string(dirname + string(theta_of_tracks->GetName())+".png").data());
  phi_of_tracks->Write();
  phi_of_tracks->Draw();
  gPad->SaveAs(string(dirname + string(phi_of_tracks->GetName())+".png").data());
  length_of_tracks->Write();
  length_of_tracks->Draw();
  gPad->SaveAs(string(dirname + string(length_of_tracks->GetName())+".png").data());
  dQds_run_view0->Write();
  dQds_run_view0->Draw();
  gPad->SaveAs(string(dirname + string(dQds_run_view0->GetName())+".png").data());
  dQds_run_view1->Write();
  dQds_run_view1->Draw();
  gPad->SaveAs(string(dirname + string(dQds_run_view1->GetName())+".png").data());
  number_of_hits_per_tracks_view0->Write();
  number_of_hits_per_tracks_view0->Draw();
  gPad->SaveAs(string(dirname + string(number_of_hits_per_tracks_view0->GetName())+".png").data());
  number_of_hits_per_tracks_view1->Write();
  number_of_hits_per_tracks_view1->Draw();
  gPad->SaveAs(string(dirname + string(number_of_hits_per_tracks_view1->GetName())+".png").data());
  number_of_tracks_per_event->Write();
  number_of_tracks_per_event->Draw();
  gPad->SaveAs(string(dirname + string(number_of_tracks_per_event->GetName())+".png").data());
  gPad->Clear();
  
  TNamed named_nb_tracks(string("nb_tracks").data(), to_string(nb_tracks).data());
  TNamed named_nb_hits_view0(string("nb_hits_view0").data(), to_string(nb_hits_view0).data());
  TNamed named_nb_hits_view1(string("nb_hits_view1").data(), to_string(nb_hits_view1).data());
  named_nb_tracks.Write();
  named_nb_hits_view0.Write();
  named_nb_hits_view1.Write();
  
  ofile.Close();
  
  return;
}


void cuts_analysis(int run = 840, string cut_type = "Ds", string version = "June", string m_dQ = "sum", string m_ds = "3D", bool save_before_cuts = true){
  method_dQ = m_dQ;
  method_ds = m_ds;
  if(!Load_Version(version)){return;}
  cut_type = set_cuts(cut_type);
  gErrorIgnoreLevel = kError;

  TTree* tree_before_cuts;
  TTree* tree_selected_tracks;
  vector<track> tracks_before_cuts;
  vector<track> selected_tracks;
  track* single_track_before_cuts = 0;
  track* single_selected_track = 0;
  int Nentries = 0;
  int item_max = 0;
  int item_min = 1;
  
  #if verbose
  cout << "Loading tree... " << endl;
  #endif
  
  string inputfile_before_cuts = SelectTrack_Output + "before_cuts";
  string inputfile_selected_tracks = SelectTrack_Output + cut_type;
  if(!ExistTest(inputfile_before_cuts)){
    cout << "ERROR: You need the tracks before cuts to analyse the cut effect" << endl;
    return;
  }
  if(!ExistTest(inputfile_selected_tracks)){
    cout << "ERROR: folder " << inputfile_selected_tracks << " not found" << endl;
    return;
  }
  inputfile_before_cuts = inputfile_before_cuts+"/tracks_"+to_string(run)+".root";
  inputfile_selected_tracks = inputfile_selected_tracks+"/tracks_"+to_string(run)+".root";
  if(!ExistTest(inputfile_before_cuts)){
    cout << "ERROR: file " << inputfile_before_cuts << " not found" << endl;
    return;
  }
  if(!ExistTest(inputfile_selected_tracks)){
    cout << "ERROR: file " << inputfile_selected_tracks << " not found" << endl;
    return;
  }
  
  cout << "Loading tracks from " << inputfile_before_cuts << endl;
  TFile ifile_before_cuts(inputfile_before_cuts.data(), "READ");
  ifile_before_cuts.GetObject("tracks", tree_before_cuts);
  item_max = tree_before_cuts->GetEntries();
  tree_before_cuts->SetBranchAddress("track",&single_track_before_cuts);
  for(int i=0; i<tree_before_cuts->GetEntries(); i++){
    tree_before_cuts->GetEntry(i);
    tracks_before_cuts.push_back(*single_track_before_cuts);
  }
  ifile_before_cuts.Close();
  
  cout << "Loading tracks from " << inputfile_selected_tracks << endl;
  TFile ifile_selected_tracks(inputfile_selected_tracks.data(), "READ");
  ifile_selected_tracks.GetObject("tracks", tree_selected_tracks);
  item_min = tree_selected_tracks->GetEntries();
  tree_selected_tracks->SetBranchAddress("track",&single_selected_track);
  for(int i=0; i<tree_selected_tracks->GetEntries(); i++){
    tree_selected_tracks->GetEntry(i);
    selected_tracks.push_back(*single_selected_track);
  }
  ifile_selected_tracks.Close();
  
  if(item_min > item_max){
    cout << "ERROR: you have more tracks after cut than before!" << endl;
    return;
  }
  if( item_min == 0 ){
    cout << "ERROR: no track found after cut" << endl;
    return;
  }
  
  #if verbose
  cout << "Analysing tracks after cut " << cut_type << "..." << endl;
  #endif
  
  string outpath = cuts_analysis_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type;
  check_and_mkdir(outpath);
  outpath = outpath + "/" + to_string(run) + "/";
  check_and_mkdir(outpath);
  analyse_cuts(selected_tracks, outpath);
  plot_selected_tracks(selected_tracks, outpath);
  plot_selected_tracks(selected_tracks, outpath, 100);

  cout << "Done." << endl;
  
  if(!save_before_cuts){return;}
  
  #if verbose
  cout << "Analysing tracks before cut..." << endl;
  #endif
  
  outpath = cuts_analysis_Output + "before_cuts";
  check_and_mkdir(outpath);
  outpath = outpath + "/" + to_string(run) + "/";
  check_and_mkdir(outpath);
  analyse_cuts(tracks_before_cuts, outpath);
  plot_selected_tracks(tracks_before_cuts, outpath);
  plot_selected_tracks(tracks_before_cuts, outpath, 100);
  
  cout << "Done." << endl;
  return;
}
