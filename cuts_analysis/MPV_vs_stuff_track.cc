#include <311Lib.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;



void MPV_vs_stuff_track(int run = 840, string cut_type = "before_cuts", string method = "length", string version = "June", string m_dQ = "sum", string m_ds = "3D" ){
  method_ds = m_ds;
  method_dQ = m_dQ;
  if(!Load_Version(version)){return;}
  cut_type = set_cuts(cut_type);
  bool MC = false;
  if(version == "MC"){
    MC = true;
    version = "June";
  }
  gStyle->SetOptFit(1111);
  gErrorIgnoreLevel = kError;
  
  #if verbose
  cout << "Doing MPV vs track " << method << " for run " << run << endl;
  #endif
  string inputfile = path_311data + "selected_tracks_" + method_ds + "_" + method_dQ + "/"+cut_type+"/tracks_"+to_string(run)+".root";
  #if verbose
  cout << "Opening file" << inputfile << " ..." << endl;
  #endif
  TFile ifile(inputfile.data(), "READ");
  if(!ifile.IsOpen()){
    cout << "ERROR: can not open input file " << inputfile << endl;
    return;
  }
  
  string outpath = MPV_vs_stuff_tracks_Output;
  check_and_mkdir(outpath);
  outpath = outpath + cut_type;
  check_and_mkdir(outpath);
  outpath = outpath +"/"+to_string(run);
  check_and_mkdir(outpath);
  outpath = outpath +"/"+method+"/";
  check_and_mkdir(outpath);
  check_and_mkdir(outpath+"dQds");
  string dirname = outpath;
  outpath = outpath + method + ".root";
  double stuff_min, dstuff, stuff_max; 
  int nstuff; 
  if(method == "length"){
    stuff_min = 0;
    dstuff = 10.;
    stuff_max = 100.; //max of 10.36 with cuts of 10deg on theta and phi
    nstuff = (stuff_max-stuff_min)/dstuff;
  }
  else{cout << "ERROR: unknown variable " << method << endl; return;}
  
  TTree *tree;
  ifile.GetObject("tracks", tree);
  track * ptrack = 0;
  vector<track> tracks;
  #if verbose
  cout << "  Loading trees... " << endl;
  #endif
  
  tree->SetBranchAddress("track",&ptrack);
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    tracks.push_back(*ptrack);
  }
  ifile.Close();
  
  vector<vector<TH1D> > dqds_per_stuff;
  map<int, vector<vector<TH1D> > > dqds_per_stuff_ByLEMs;
  vector<double> stuffs;
  vector<TGraphErrors> mpv_vs_stuff;
  map<int, vector<TGraphErrors> > mpv_vs_stuff_ByLEMs;
  TGraph expected_MPV;
  expected_MPV.SetName("expected_MPV");
  expected_MPV.SetLineColor(kBlue);
  init_graph_purity_and_charging_up(mpv_vs_stuff, mpv_vs_stuff_ByLEMs, string("mpv_vs_"+method));
  for(int i = 0; i < nstuff; i++){
    stuffs.push_back(stuff_min+i*dstuff);
    string str_stuff = to_string_with_precision( stuff_min+i*dstuff, 4 );
    if( str_stuff.find(".") != string::npos ){str_stuff.replace(str_stuff.find("."),1,"_");}
    
    dqds_per_stuff.push_back({TH1D(string("dqds_view_0_at_"+method+"_"+str_stuff).data(),string("dqds_view_0_at_"+method+"_"+str_stuff).data(),100,0,50), TH1D(string("dqds_view_1_at_"+method+"_"+str_stuff).data(),string("dqds_view_1_at_"+method+"_"+str_stuff).data(),100,0,50)});
    
    for( auto lem : lems ){
      dqds_per_stuff_ByLEMs[lem].push_back({TH1D(string("dqds_view_0_LEM_"+to_string(lem)+"_at_"+method+"_"+str_stuff).data(),string("dqds_view_0_LEM_"+to_string(lem)+"_at_"+method+"_"+str_stuff).data(),100,0,50), TH1D(string("dqds_view_1_LEM_"+to_string(lem)+"_at_"+method+"_"+str_stuff).data(),string("dqds_view_1_LEM_"+to_string(lem)+"_at_"+method+"_"+str_stuff).data(),100,0,50)});
    }
  }
  for( auto t : tracks){
    for( auto h : t.hits_trk ){
      if( !isGood_lem(h.lem)){continue;}
      double dq;
      double ds;
      if(method_ds == "3D"){ds = h.dx_3D;}
      else{ds = h.dx_local;}
      if(method_dQ == "sum"){dq = h.dq_sum;}
      else{dq = h.dq_integral;}
      double dqds = h.gain_density_correction_factor * dq/ds;
      //cut on dqdx
      if( dqds <= 0 || dqds > 50 || h.sp_x < tpc_boundaries[0] || h.sp_x >= tpc_boundaries[1] || h.sp_y < tpc_boundaries[2] || h.sp_y >= tpc_boundaries[3] || h.sp_z < tpc_boundaries[4] || h.sp_z >= tpc_boundaries[5] ){ continue; }
      int istuff = -1;
      if(method == "length" && !MC){istuff = (int)((t.length-stuff_min)/dstuff);}
      if(istuff < 0){continue;}
      if(istuff > nstuff-1){continue;}
      dqds_per_stuff[istuff][h.view].Fill(dqds);
      dqds_per_stuff_ByLEMs[h.lem][istuff][h.view].Fill(dqds);
    }
  }
  
  /////////////////////////////////////////////////////////////////
  
  TFile ofile(outpath.data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << "ERROR: can not open output file " << outpath << endl;
    return;
  }
  
  for(int i = 0; i < nstuff; i++){
    if(method == "ds"){
      expected_MPV.SetPoint(expected_MPV.GetN(),stuffs[i]+0.5*dstuff,MPVs_from_theory(stuffs[i]+0.5*dstuff)/2.);
    }
    vector<double> f0 = {-1,-1};
    vector<double> f1 = {-1,-1};
    f0 = fit_dQds(&dqds_per_stuff[i][0], false, 200, 0.05, 4, &mpv_vs_stuff[0], stuffs[i]+0.5*dstuff, 0.5*dstuff);
    dqds_per_stuff[i][0].Write();
    dqds_per_stuff[i][0].Draw();
    gPad->SaveAs(string(dirname + "dQds/" + string(dqds_per_stuff[i][0].GetName())+".png").data());
    f1 = fit_dQds(&dqds_per_stuff[i][1], false, 200, 0.05, 4, &mpv_vs_stuff[1], stuffs[i]+0.5*dstuff, 0.5*dstuff);
    dqds_per_stuff[i][1].Write();
    dqds_per_stuff[i][1].Draw();
    gPad->SaveAs(string(dirname + "dQds/" + string(dqds_per_stuff[i][1].GetName())+".png").data());
    for( auto lem : lems ){
      f0 = fit_dQds(&dqds_per_stuff_ByLEMs[lem][i][0], false, 200, 0.05, 4, &mpv_vs_stuff_ByLEMs[lem][0], stuffs[i]+0.5*dstuff, 0.5*dstuff);
      f1 = fit_dQds(&dqds_per_stuff_ByLEMs[lem][i][1], false, 200, 0.05, 4, &mpv_vs_stuff_ByLEMs[lem][1], stuffs[i]+0.5*dstuff, 0.5*dstuff);
      dqds_per_stuff_ByLEMs[lem][i][0].Write();
      dqds_per_stuff_ByLEMs[lem][i][1].Write();
      if( f0[0] > 0 && f0[0] > 0){
        if(mpv_cosmics > 0){
          f0[0] = f0[0]/mpv_cosmics;
          f0[1] = f0[1]/mpv_cosmics;
          f1[0] = f1[0]/mpv_cosmics;
          f1[1] = f1[1]/mpv_cosmics;
        }
        mpv_vs_stuff_ByLEMs[lem][2].SetPoint(mpv_vs_stuff_ByLEMs[lem][2].GetN(),stuffs[i]+0.5*dstuff,f0[0]+f1[0]);
        mpv_vs_stuff_ByLEMs[lem][2].SetPointError(mpv_vs_stuff_ByLEMs[lem][2].GetN()-1,0.5*dstuff,TMath::Sqrt(f0[1]*f0[1]+f1[1]*f1[1]));
      }
    }
  }
  
  
  mpv_vs_stuff[0].SetName(string("mpv_vs_"+method+"view_0").data());
  mpv_vs_stuff[0].Draw("AP");
  mpv_vs_stuff[0].SetMinimum(0);
  mpv_vs_stuff[0].SetMaximum(15);
  mpv_vs_stuff[0].GetXaxis()->SetRangeUser(stuff_min,stuff_max);
  if(method == "ds"){
    expected_MPV.Draw("SAME");
  }
  ofile.cd();
  mpv_vs_stuff[0].Write();
  gPad->SaveAs(string(dirname + string(mpv_vs_stuff[0].GetName())+".png").data());

  mpv_vs_stuff[1].SetName(string("mpv_vs_"+method+"view_1").data());
  mpv_vs_stuff[1].Draw("AP");
  mpv_vs_stuff[1].SetMinimum(0);
  mpv_vs_stuff[1].SetMaximum(15);
  mpv_vs_stuff[1].GetXaxis()->SetRangeUser(stuff_min,stuff_max);
  if(method == "ds"){
    expected_MPV.Draw("SAME");
  }
  ofile.cd();
  mpv_vs_stuff[1].Write();
  gPad->SaveAs(string(dirname + string(mpv_vs_stuff[1].GetName())+".png").data());

  mpv_vs_stuff[2].SetName(string("summed_views_mpv_vs_"+method).data());
  mpv_vs_stuff[2].Draw("AP");
  mpv_vs_stuff[2].SetMinimum(-5);
  mpv_vs_stuff[2].SetMaximum(5);
  mpv_vs_stuff[2].GetXaxis()->SetRangeUser(stuff_min,stuff_max);
  if(method == "ds"){
    expected_MPV.Draw("SAME");
  }
  ofile.cd();
  mpv_vs_stuff[2].Write();
  gPad->SaveAs(string(dirname + string(mpv_vs_stuff[2].GetName())+".png").data());
  
  gPad->Clear();
  
  for( auto lem : lems ){
    mpv_vs_stuff_ByLEMs[lem][0].SetName(string("mpv_vs_"+method+"_lem_"+to_string(lem)+"view_0").data());
    mpv_vs_stuff_ByLEMs[lem][0].Draw("AP");
    mpv_vs_stuff_ByLEMs[lem][0].SetMinimum(0);
    mpv_vs_stuff_ByLEMs[lem][0].SetMaximum(15);
    mpv_vs_stuff_ByLEMs[lem][0].GetXaxis()->SetRangeUser(stuff_min,stuff_max);
    if(method == "ds"){
      expected_MPV.Draw("SAME");
    }
    ofile.cd();
    mpv_vs_stuff_ByLEMs[lem][0].Write();
    gPad->SaveAs(string(dirname + string(mpv_vs_stuff_ByLEMs[lem][0].GetName())+".png").data());

    mpv_vs_stuff_ByLEMs[lem][1].SetName(string("mpv_vs_"+method+"_lem_"+to_string(lem)+"view_1").data());
    mpv_vs_stuff_ByLEMs[lem][1].Draw("AP");
    mpv_vs_stuff_ByLEMs[lem][1].SetMinimum(0);
    mpv_vs_stuff_ByLEMs[lem][1].SetMaximum(15);
    mpv_vs_stuff_ByLEMs[lem][1].GetXaxis()->SetRangeUser(stuff_min,stuff_max);
    if(method == "ds"){
      expected_MPV.Draw("SAME");
    }
    ofile.cd();
    mpv_vs_stuff_ByLEMs[lem][1].Write();
    gPad->SaveAs(string(dirname + string(mpv_vs_stuff_ByLEMs[lem][1].GetName())+".png").data());

    mpv_vs_stuff_ByLEMs[lem][2].SetName(string("summed_views_mpv_vs_"+method+"_lem_"+to_string(lem)).data());
    mpv_vs_stuff_ByLEMs[lem][2].Draw("AP");
    mpv_vs_stuff_ByLEMs[lem][2].SetMinimum(-5);
    mpv_vs_stuff_ByLEMs[lem][2].SetMaximum(5);
    mpv_vs_stuff_ByLEMs[lem][2].GetXaxis()->SetRangeUser(stuff_min,stuff_max);
    if(method == "ds"){
      expected_MPV.Draw("SAME");
    }
    ofile.cd();
    mpv_vs_stuff_ByLEMs[lem][2].Write();
    gPad->SaveAs(string(dirname + string(mpv_vs_stuff_ByLEMs[lem][2].GetName())+".png").data());
    
    gPad->Clear();
  }
  ofile.Close();

  return;
}
