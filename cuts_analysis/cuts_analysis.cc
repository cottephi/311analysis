#include <311Lib.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}

void analyse_cuts(vector<track> tracks, string cut_type, double dqds_cut = 0){

  int run = tracks[0].run;
  
  map<int,int> events_and_number_of_tracks;
  int nb_tracks = 0;
  int nb_hits_view0 = 0;
  int nb_hits_view1 = 0;
  
  int max_hits_per_event = 0;
  int max_tracks_per_event = 0;
  
  TH1D dQds_run_view0("dQds_run_view0","dQds_run_view0",100,0,50);
  TH1D dQds_run_view1("dQds_run_view1","dQds_run_view1",100,0,50);
  TH1D dQds_run_view0_20("dQds_run_view0_20","dQds_run_view0_20",100,0,20);
  TH1D dQds_run_view1_20("dQds_run_view1_20","dQds_run_view1_20",100,0,20);
  TH1D ds_of_hits_view0("ds_of_hits_view0","ds_of_hits_view0",100,0,6);
  TH1D ds_of_hits_view1("ds_of_hits_view1","ds_of_hits_view1",100,0,6);
  TH1D GoF_of_hits_view0("GoF_of_hits_view0","GoF_of_hits_view0",100,0,10);
  TH1D GoF_of_hits_view1("GoF_of_hits_view1","GoF_of_hits_view1",100,0,10);
  
  TH1D dQds_theta_0_180_to_160("dQds_theta_0_180_to_160","dQds_theta_0_180_to_160",100,0,50);
  TH1D dQds_theta_0_180_to_140("dQds_theta_0_180_to_140","dQds_theta_0_180_to_140",100,0,50);
  TH1D dQds_theta_0_180_to_120("dQds_theta_0_180_to_120","dQds_theta_0_180_to_120",100,0,50);
  TH1D dQds_theta_0_180_to_100("dQds_theta_0_180_to_100","dQds_theta_0_180_to_100",100,0,50);
  TH1D dQds_theta_0_180_to_90("dQds_theta_0_180_to_90","dQds_theta_0_180_to_90",100,0,50);
  
  TH1D dQds_theta_1_180_to_160("dQds_theta_1_180_to_160","dQds_theta_1_180_to_160",100,0,50);
  TH1D dQds_theta_1_180_to_140("dQds_theta_1_180_to_140","dQds_theta_1_180_to_140",100,0,50);
  TH1D dQds_theta_1_180_to_120("dQds_theta_1_180_to_120","dQds_theta_1_180_to_120",100,0,50);
  TH1D dQds_theta_1_180_to_100("dQds_theta_1_180_to_100","dQds_theta_1_180_to_100",100,0,50);
  TH1D dQds_theta_1_180_to_90("dQds_theta_1_180_to_90","dQds_theta_1_180_to_90",100,0,50);
  
  TH2D dQds_vs_hit_GoF_view0("dQds_vs_hit_GoF_view0","dQds_vs_hit_GoF_view0",100,0,10,100,0,50);
  TH2D dQds_vs_hit_GoF_view1("dQds_vs_hit_GoF_view1","dQds_vs_hit_GoF_view1",100,0,10,100,0,50);
  TH2D dQds_vs_hit_ds_view0("dQds_vs_hit_ds_view0","dQds_vs_hit_ds_view0",100,0,6,100,0,50);
  TH2D dQds_vs_hit_ds_view1("dQds_vs_hit_ds_view1","dQds_vs_hit_ds_view1",100,0,6,100,0,50);
  TH2D dQds_vs_hit_dQ_view0("dQds_vs_hit_dQ_view0","dQds_vs_hit_dQ_view0",100,0,200,100,0,50);
  TH2D dQds_vs_hit_dQ_view1("dQds_vs_hit_dQ_view1","dQds_vs_hit_dQ_view1",100,0,200,100,0,50);
  TH2D ds_local_vs_ds_phil_view0("ds_local_vs_ds_phil_view0","ds_local_vs_ds_phil_view0",100,0,50,100,0,50);
  TH2D ds_local_vs_ds_phil_view1("ds_local_vs_ds_phil_view1","ds_local_vs_ds_phil_view1",100,0,50,100,0,50);
  TH2D sumQ0_vs_sumQ1_1000("sumQ0_vs_sumQ1_1000","sumQ0_vs_sumQ1_1000",100,0,1000,100,0,1000);
  TH2D sumQ0_vs_sumQ1_2000("sumQ0_vs_sumQ1_2000","sumQ0_vs_sumQ1_2000",100,0,2000,100,0,2000);
  TH2D sumQ0_vs_sumQ1("sumQ0_vs_sumQ1","sumQ0_vs_sumQ1",100,0,10000,100,0,10000);
  
  TH1D length_of_tracks("length_of_tracks","length_of_tracks",100,0,332);
  TH1D theta_of_tracks("theta_of_tracks","theta_of_tracks",100,90,180);
  TH1D phi_of_tracks("phi_of_tracks","phi_of_tracks",100,-180,180);
  TH1D number_of_hits_per_tracks_view0("number_of_hits_per_tracks_view0","tracks_and_number_of_hits_per_tracks_view0",100,0,1120);
  TH1D number_of_hits_per_tracks_view1("number_of_hits_per_tracks_view1","tracks_and_number_of_hits_per_tracks_view1",100,0,1120);
  TH1D number_of_tracks_per_event("number_of_tracks_per_event","number_of_tracks_per_event",71,0,70);
  TH1D ratioQ1overQ0("ratioQ1overQ0","ratioQ1overQ0",100,0,2);
  
  TH2D dQds_vs_track_length_view0("dQds_vs_track_length_view0","dQds_vs_track_length_view0",100,0,300,100,0,50);
  TH2D dQds_vs_track_length_view1("dQds_vs_track_length_view1","dQds_vs_track_length_view1",100,0,300,100,0,50);
  TH2D dQds_vs_ntracksinevent_view0("dQds_vs_ntracksinevent_view0","dQds_vs_ntracksinevent_view0",71,1,70,100,0,50);
  TH2D dQds_vs_ntracksinevent_view1("dQds_vs_ntracksinevent_view1","dQds_vs_ntracksinevent_view1",71,1,70,100,0,50);
  TH2D dQds_vs_track_theta_view0_0("dQds_vs_track_theta_view0_0","dQds_vs_track_theta_view0",100,90,105,100,0,50);
  TH2D dQds_vs_track_theta_view1_0("dQds_vs_track_theta_view1_0","dQds_vs_track_theta_view1",100,90,105,100,0,50);
  TH2D dQds_vs_track_theta_view0_1("dQds_vs_track_theta_view0_1","dQds_vs_track_theta_view0",100,105,90+160,100,0,50);
  TH2D dQds_vs_track_theta_view1_1("dQds_vs_track_theta_view1_1","dQds_vs_track_theta_view1",100,105,160,100,0,50);
  TH2D dQds_vs_track_theta_view0_2("dQds_vs_track_theta_view0_2","dQds_vs_track_theta_view0",100,160,180,100,0,50);
  TH2D dQds_vs_track_theta_view1_2("dQds_vs_track_theta_view1_2","dQds_vs_track_theta_view1",100,160,180,100,0,50);
  TH2D dQds_vs_track_phi_view0("dQds_vs_track_phi_view0","dQds_vs_track_phi_view0",100,-180,180,100,0,50);
  TH2D dQds_vs_track_phi_view1("dQds_vs_track_phi_view1","dQds_vs_track_phi_view1",100,-180,180,100,0,50);
  TH2D dQds_vs_track_ds_view0("dQds_vs_track_ds_view0","dQds_vs_track_ds_view0",100,0,1,100,0,50);
  TH2D dQds_vs_track_ds_view1("dQds_vs_track_ds_view1","dQds_vs_track_ds_view1",100,0,1,100,0,50);
  TH2D dQds_vs_number_of_hits_per_tracks_view0("dQds_vs_number_of_hits_per_tracks_view0","dQds_vs_number_of_hits_per_tracks_view0",100,1,100,100,0,50);
  TH2D dQds_vs_number_of_hits_per_tracks_view1("dQds_vs_number_of_hits_per_tracks_view1","dQds_vs_number_of_hits_per_tracks_view1",100,1,100,100,0,50);

  TH2D xz("xz","xz",300, 0, 300,100,-50,50);
  TH2D xy("xy","xy",100,-50,50, 100,-50,50);
  TH2D yz("yz","yz",300, 0, 300,100,-50,50);
  TH3D xyz("xyz","xyz",300,0,300,100,-50,50,100,-50,50);

  TH1D ds_of_tracks_view0("ds_of_tracks_view0","ds_of_tracks_view0",100,0,10);
  TH1D ds_of_tracks_view1("ds_of_tracks_view1","ds_of_tracks_view1",100,0,10);
  
  double ds;
  double dq;
  double dqds;
  
  for( auto t : tracks){
//    #if verbose
//    cout << "Analysing track " << t.event << " " << t.id << endl;
//    #endif
    nb_tracks++;
    int hits_track_view0 = 0;
    int hits_track_view1 = 0;
    double sumQ0 = 0;
    double sumQ1 = 0;
    
    for(auto h : t.hits_trk){
      if(h.view == 0){hits_track_view0++;}
      if(h.view == 1){hits_track_view1++;}
    }

    for(auto h : t.hits_trk){

      if(method_ds == "3D"){ds = h.dx_3D;}
      else if(method_ds == "local"){ds = h.dx_local;}
      else{ds = h.dx_phil;}
      if(method_dQ == "sum"){dq = h.dq_sum;}
      else{dq = h.dq_integral;}
      dqds = dq/ds;
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
      if(h.view == 0){
        sumQ0 += dq;
        dQds_run_view0.Fill(dqds);
        dQds_run_view0_20.Fill(dqds);
        ds_of_hits_view0.Fill(ds);
        GoF_of_hits_view0.Fill(h.GoF);
        nb_hits_view0++;
        dQds_vs_track_length_view0.Fill(t.length,dqds);
        dQds_vs_ntracksinevent_view0.Fill(t.event_ntracks,dqds);
        dQds_vs_track_theta_view0_0.Fill(t.theta,dqds);
        dQds_vs_track_theta_view0_1.Fill(t.theta,dqds);
        dQds_vs_track_theta_view0_2.Fill(t.theta,dqds);
        if(t.theta > 90){
          dQds_theta_0_180_to_90.Fill(dqds);
        }
        if(t.theta > 100){
          dQds_theta_0_180_to_100.Fill(dqds);
        }
        if(t.theta > 120){
          dQds_theta_0_180_to_120.Fill(dqds);
        }
        if(t.theta > 140){
          dQds_theta_0_180_to_140.Fill(dqds);
        }
        if(t.theta > 160){
          dQds_theta_0_180_to_160.Fill(dqds);
        }
        dQds_vs_track_phi_view0.Fill(t.phi,dqds);
        dQds_vs_track_ds_view0.Fill(get_dss(t)[0],dqds);
        dQds_vs_number_of_hits_per_tracks_view0.Fill(hits_track_view0,dqds);
        dQds_vs_hit_GoF_view0.Fill(h.GoF,dqds);
        dQds_vs_hit_ds_view0.Fill(ds,dqds);
        dQds_vs_hit_dQ_view0.Fill(dq,dqds);
        ds_local_vs_ds_phil_view0.Fill(h.dx_local,h.dx_phil);
      }
      if(h.view == 1){
        sumQ1 += dq;
        dQds_run_view1.Fill(dqds);
        dQds_run_view1_20.Fill(dqds);
        ds_of_hits_view1.Fill(ds);
        GoF_of_hits_view1.Fill(h.GoF);
        hits_track_view1++;
        nb_hits_view1++;
        dQds_vs_track_length_view1.Fill(t.length,dqds);
        dQds_vs_ntracksinevent_view1.Fill(t.event_ntracks,dqds);
        dQds_vs_track_theta_view1_0.Fill(t.theta,dqds);
        dQds_vs_track_theta_view1_1.Fill(t.theta,dqds);
        dQds_vs_track_theta_view1_2.Fill(t.theta,dqds);
        if(t.theta > 90){
          dQds_theta_1_180_to_90.Fill(dqds);
        }
        if(t.theta > 100){
          dQds_theta_1_180_to_100.Fill(dqds);
        }
        if(t.theta > 120){
          dQds_theta_1_180_to_120.Fill(dqds);
        }
        if(t.theta > 140){
          dQds_theta_1_180_to_140.Fill(dqds);
        }
        if(t.theta > 160){
          dQds_theta_1_180_to_160.Fill(dqds);
        }
        dQds_vs_track_phi_view1.Fill(t.phi,dqds);
        dQds_vs_track_ds_view1.Fill(get_dss(t)[0],dqds);
        dQds_vs_number_of_hits_per_tracks_view1.Fill(hits_track_view1,dqds);
        dQds_vs_hit_GoF_view1.Fill(h.GoF,dqds);
        dQds_vs_hit_ds_view1.Fill(ds,dqds);
        dQds_vs_hit_dQ_view1.Fill(dq,dqds);
        ds_local_vs_ds_phil_view1.Fill(h.dx_local,h.dx_phil);
      }
    }
    number_of_tracks_per_event.Fill(t.event_ntracks);
    xz.Fill(t.start_z,t.start_x);
    xz.Fill(t.end_z,t.end_x);
    xy.Fill(t.start_y,t.start_x);
    xy.Fill(t.end_y,t.end_x);
    yz.Fill(t.start_z,t.start_y);
    yz.Fill(t.end_z,t.end_y);
    xyz.Fill(t.start_z,t.start_y,t.start_x);
    xyz.Fill(t.end_z,t.end_y,t.end_x);
    number_of_hits_per_tracks_view0.Fill(hits_track_view0);
    number_of_hits_per_tracks_view1.Fill(hits_track_view1);
    length_of_tracks.Fill(t.length);
    theta_of_tracks.Fill(t.theta);
    phi_of_tracks.Fill(t.phi);
    ds_of_tracks_view0.Fill(get_dss(t)[0]);
    ds_of_tracks_view1.Fill(get_dss(t)[1]);
    sumQ0_vs_sumQ1_1000.Fill(sumQ0,sumQ1);
    sumQ0_vs_sumQ1_2000.Fill(sumQ0,sumQ1);
    sumQ0_vs_sumQ1.Fill(sumQ0,sumQ1);
    ratioQ1overQ0.Fill(sumQ1/sumQ0);
  }
  
  
  ///////////////////////////////////////save///////////////////////////////////////////
  
  string outpath = cuts_analysis_Output + cut_type + "/" + to_string(run) + "/";
  check_and_mkdir(outpath);
  if(dqds_cut != 0){
    if(dqds_cut < 0){
      outpath = outpath + "dQds_inf_to_" + to_string_with_precision(dqds_cut,2);
    }
    else{
      outpath = outpath + "dQds_sup_to_" + to_string_with_precision(dqds_cut,2);
    }
  }
  string outfile = outpath + "analysis.root";
  TFile ofile(outfile.data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << "ERROR: can not open " << outpath << endl;
    return;
  }
  #if verbose
  cout << " TFile " << outfile << " has been created " << endl;
  #endif
  
  TCanvas mycan("mycan","mycan",1000,618);
  dQds_run_view0.Write();
  dQds_run_view0.Draw();
  mycan.Print(string(outpath+"analysis.pdf(").data(),"pdf");
  dQds_run_view1.Write();
  dQds_run_view1.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_run_view0_20.Write();
  dQds_run_view0_20.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_run_view1_20.Write();
  dQds_run_view1_20.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ds_of_tracks_view0.Write();
  ds_of_tracks_view0.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ds_of_tracks_view1.Write();
  ds_of_tracks_view1.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  xyz.Write();
  xyz.Draw("COL");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  xz.Write();
  xz.Draw("COL");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  xy.Write();
  xy.Draw("COL");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  yz.Write();
  yz.Draw("COL");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ds_of_hits_view0.Write();
  ds_of_hits_view0.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ds_of_hits_view1.Write();
  ds_of_hits_view1.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  GoF_of_hits_view0.Write();
  GoF_of_hits_view0.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  GoF_of_hits_view1.Write();
  GoF_of_hits_view1.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  theta_of_tracks.Write();
  theta_of_tracks.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  phi_of_tracks.Write();
  phi_of_tracks.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  length_of_tracks.Write();
  length_of_tracks.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  number_of_hits_per_tracks_view0.Write();
  number_of_hits_per_tracks_view0.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  number_of_hits_per_tracks_view1.Write();
  number_of_hits_per_tracks_view1.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  number_of_tracks_per_event.Write();
  number_of_tracks_per_event.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_ntracksinevent_view0.Write();
  dQds_vs_ntracksinevent_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_ntracksinevent_view1.Write();
  dQds_vs_ntracksinevent_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_length_view0.Write();
  dQds_vs_track_length_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_length_view1.Write();
  dQds_vs_track_length_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_theta_view0_0.Write();
  dQds_vs_track_theta_view0_0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_theta_view0_1.Write();
  dQds_vs_track_theta_view0_1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_theta_view0_2.Write();
  dQds_vs_track_theta_view0_2.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_0_180_to_90.Write();
  dQds_theta_0_180_to_90.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_0_180_to_100.Write();
  dQds_theta_0_180_to_100.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_0_180_to_120.Write();
  dQds_theta_0_180_to_120.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_0_180_to_140.Write();
  dQds_theta_0_180_to_140.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_0_180_to_160.Write();
  dQds_theta_0_180_to_160.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_theta_view1_0.Write();
  dQds_vs_track_theta_view1_0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_theta_view1_1.Write();
  dQds_vs_track_theta_view1_1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_theta_view1_2.Write();
  dQds_vs_track_theta_view1_2.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_1_180_to_90.Write();
  dQds_theta_1_180_to_90.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_1_180_to_100.Write();
  dQds_theta_1_180_to_100.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_1_180_to_120.Write();
  dQds_theta_1_180_to_120.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_1_180_to_140.Write();
  dQds_theta_1_180_to_140.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_theta_1_180_to_160.Write();
  dQds_theta_1_180_to_160.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_phi_view0.Write();
  dQds_vs_track_phi_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_phi_view1.Write();
  dQds_vs_track_phi_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_ds_view0.Write();
  dQds_vs_track_ds_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_track_ds_view1.Write();
  dQds_vs_track_ds_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_number_of_hits_per_tracks_view0.Write();
  dQds_vs_number_of_hits_per_tracks_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_number_of_hits_per_tracks_view1.Write();
  dQds_vs_number_of_hits_per_tracks_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_hit_GoF_view0.Write();
  dQds_vs_hit_GoF_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_hit_GoF_view1.Write();
  dQds_vs_hit_GoF_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_hit_ds_view0.Write();
  dQds_vs_hit_ds_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_hit_ds_view1.Write();
  dQds_vs_hit_ds_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_hit_dQ_view0.Write();
  dQds_vs_hit_dQ_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  dQds_vs_hit_dQ_view1.Write();
  dQds_vs_hit_dQ_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ds_local_vs_ds_phil_view0.Write();
  ds_local_vs_ds_phil_view0.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ds_local_vs_ds_phil_view1.Write();
  ds_local_vs_ds_phil_view1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  ratioQ1overQ0.Write();
  ratioQ1overQ0.Draw();
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  sumQ0_vs_sumQ1_1000.Write();
  sumQ0_vs_sumQ1_1000.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  sumQ0_vs_sumQ1_2000.Write();
  sumQ0_vs_sumQ1_2000.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf").data(),"pdf");
  sumQ0_vs_sumQ1.Write();
  sumQ0_vs_sumQ1.Draw("COLZ");
  mycan.Print(string(outpath+"analysis.pdf)").data(),"pdf");
  
  TNamed named_nb_tracks(string("nb_tracks").data(), to_string(nb_tracks).data());
  TNamed named_nb_hits_view0(string("nb_hits_view0").data(), to_string(nb_hits_view0).data());
  TNamed named_nb_hits_view1(string("nb_hits_view1").data(), to_string(nb_hits_view1).data());
  named_nb_tracks.Write();
  named_nb_hits_view0.Write();
  named_nb_hits_view1.Write();
  
  ofile.Close();
  
  return;
}

// /*1165, 1166, 1167, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1180, 1181, */1182, 1183, 1187, 1188, 1189, 1190, 1191, /*1192, 1193,*/ 1194, 1195, 1196, 1197, /*1198, 1199,*/ 840
void cuts_analysis(vector<int> run_list={1197}, string cut_type = "philippe_theta", string v = "July", string m_dQ = "sum", string m_ds = "local"){
  method_dQ = m_dQ;
  method_ds = m_ds;
  if(!Load_Version(v)){return;}
  cut_type = set_cuts(cut_type);
  gErrorIgnoreLevel = kError;
  
  string path = path_wa105_311data;
  if (run_list.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path_wa105_311data << "*..." << endl;
    #endif
    string wildcard_path = path_wa105_311data + "*";
    for( auto irun : glob(wildcard_path) ){
      run_list.push_back(atoi(irun.data()));
    }
  }
  
  for(auto run : run_list){

    TTree* tree;
    vector<track> selected_tracks;
    track* single_selected_track = 0;
    
    #if verbose
    cout << "Loading tree... " << endl;
    #endif
    
    string inputfile = SelectTrack_Output + cut_type;
    if(!ExistTest(inputfile)){
      cout << "ERROR: folder " << inputfile << " not found" << endl;
      return;
    }
    inputfile = inputfile+"/tracks_"+to_string(run)+".root";
    if(!ExistTest(inputfile)){
      cout << "ERROR: file " << inputfile << " not found" << endl;
      return;
    }
    
    cout << "Loading tracks from " << inputfile << endl;
    TFile ifile(inputfile.data(), "READ");
    ifile.GetObject("tracks", tree);
    if( tree->GetEntries() == 0 ){
      cout << "ERROR: no tracks found" << endl;
      return;
    }
    tree->SetBranchAddress("track",&single_selected_track);
    for(int i=0; i<tree->GetEntries(); i++){
      tree->GetEntry(i);
      selected_tracks.push_back(*single_selected_track);
    }
    ifile.Close();
    
    
    #if verbose
    cout << "Analysing tracks for cut " << cut_type << " for run " << run << "..." << endl;
    #endif
    
    analyse_cuts(selected_tracks, cut_type);
  //  analyse_cuts(selected_tracks, cut_type, 10);
  //  analyse_cuts(selected_tracks, cut_type, -10);
    plot_tracks(selected_tracks, cut_type,100);
  //  plot_tracks(selected_tracks, cut_type, 100, 10);
  //  plot_tracks(selected_tracks, cut_type, 100, -10);
  }

  cout << "Done." << endl;
  return;
}
