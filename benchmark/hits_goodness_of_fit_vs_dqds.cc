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

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3){
  ostringstream out;
  out << setprecision(n) << a_value;
  return out.str();
}


void hits_goodness_of_fit_vs_dqds(int run = 840, string version = "June"){

  if(!Load_Version(version)){return;}

  gErrorIgnoreLevel = kError;
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.35);                
  // Set y-position (fraction of pad size)
  gStyle->SetStatX(0.85);                
  // Set x-position (fraction of pad size)
  gStyle->SetStatW(0.2);                
  // Set width of stat-box (fraction of pad size)
  gStyle->SetStatH(0.15);                
  // Set height of stat-box (fraction of pad size)
  string path = SelectTrack_Output;
  if(!load_run_lists()){return;}
  vector<string> run_list;
  load_cosmics();
  int min_number_of_hits = 0;
  int time = 1;

  //****************************************************************************
  //variables definition here
  
//Generated with python macro (all runs)
  
  
  TChain chain("analysistree/anatree");
  string foldername = path_wa105_311data+to_string(run)+"/";
  int n_subruns = glob(string(foldername+"*").data()).size();
  if( n_subruns == 0){
    cout << foldername << " not found or empty" << endl;
    return;
  }
  #if verbose
  cout << "Processing " << n_subruns << " subruns in run " << run << endl;
  #endif

  string file = "";
  string filename = "";
  for(int i=0; i< n_subruns; i++){
    
    file=to_string(run)+"-"+to_string(i);
    if(file == "1009-21"){continue;}
    if(path_311data.find("Feb") != string::npos){
      file = file + "-Parser.root";
    }
    else if(path_311data.find("June") != string::npos){
      file = file + "-RecoFast-Parser.root";
    }
    else{
      cout << "ERROR: unknown reconstruction version" << endl;
      return;
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
      continue;
    }
  }
  #if verbose
  cout << "Files added. " << endl;
  #endif
  
  #if verbose
  cout << "Plotting hits goodness of fit for run " << run << endl;
  #endif
  
  //multigraph for mpv, containing view 0 and view 1
  vector<TH2D> HGoF;
  HGoF.push_back(TH2D("HGoF_0","HGoF_0",100,0,50,100,-1,10));
  HGoF.push_back(TH2D("HGoF_1","HGoF_1",100,0,50,100,-1,10));
  HGoF.push_back(TH2D("HGoF_summed_views","HGoF_summed_views",100,0,50,100,-1,10));
  map<int, vector<TH2D>  > HGoF_ByLEMs;
  vector<TH2D> HGoF_integral_local;
  HGoF_integral_local.push_back(TH2D("HGoF_integral_local_integral_local_0","HGoF_integral_local_0",100,0,50,100,-1,10));
  HGoF_integral_local.push_back(TH2D("HGoF_integral_local_1","HGoF_integral_local_1",100,0,50,100,-1,10));
  HGoF_integral_local.push_back(TH2D("HGoF_integral_local_summed_views","HGoF_integral_local_summed_views",100,0,50,100,-1,10));
  map<int, vector<TH2D>  > HGoF_integral_local_ByLEMs;
  #if verbose 
  cout << "Initialising histos..." << endl;
  #endif
  for(auto lem : lems){
    HGoF_ByLEMs[lem].push_back(TH2D(string("HGoF_0_LEM_"+to_string(lem)).data(),string("HGoF_0_LEM_"+to_string(lem)).data(),100,0,50,100,-1,10));
    HGoF_ByLEMs[lem].push_back(TH2D(string("HGoF_1_LEM_"+to_string(lem)).data(),string("HGoF_1_LEM_"+to_string(lem)).data(),100,0,50,100,-1,10));
    HGoF_ByLEMs[lem].push_back(TH2D(string("HGoF_summed_views_LEM_"+to_string(lem)).data(),string("HGoF_summed_views_LEM_"+to_string(lem)).data(),100,0,50,100,-1,10));
    
    HGoF_integral_local_ByLEMs[lem].push_back(TH2D(string("HGoF_integral_local_0_LEM_"+to_string(lem)).data(),string("HGoF_integral_local_0_LEM_"+to_string(lem)).data(),100,0,50,100,-1,10));
    HGoF_integral_local_ByLEMs[lem].push_back(TH2D(string("HGoF_integral_local_1_LEM_"+to_string(lem)).data(),string("HGoF_integral_local_1_LEM_"+to_string(lem)).data(),100,0,50,100,-1,10));
    HGoF_integral_local_ByLEMs[lem].push_back(TH2D(string("HGoF_integral_local_summed_views_LEM_"+to_string(lem)).data(),string("HGoF_integral_local_summed_views_LEM_"+to_string(lem)).data(),100,0,50,100,-1,10));
  }
  const int NMaxHitsPerEvent=100000;
  const int NMaxTracksPerEvent=1000;
  int NEntries = (int)chain.GetEntries();
  int tNumberOfHits;
  short tNumberOfTracks;
  float tHit_GoodnessOfFit[NMaxHitsPerEvent];
  short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2];
  double tTrack_Hit_Y[NMaxHitsPerEvent];
  double tTrack_Hit_Z[NMaxHitsPerEvent];
  float tTrack_dx_LocalTrackDirection[NMaxHitsPerEvent];
  float tTrack_dx_3DPosition[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
  float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];
  
  chain.SetBranchAddress("NumberOfHits",&tNumberOfHits);
  chain.SetBranchAddress("Hit_GoodnessOfFit",&tHit_GoodnessOfFit);
  if(version == "June"){
    chain.SetBranchAddress("NumberOfTracks",&tNumberOfTracks);  
    chain.SetBranchAddress("Track_NumberOfHitsPerView",&tTrack_NumberOfHitsPerView);
    chain.SetBranchAddress("Track_Hit_Y", &tTrack_Hit_Y);
    chain.SetBranchAddress("Track_Hit_Z", &tTrack_Hit_Z);
  chain.SetBranchAddress("Track_Hit_ds_LocalTrackDirection", &tTrack_dx_LocalTrackDirection);
  chain.SetBranchAddress("Track_Hit_ds_3DPosition", &tTrack_dx_3DPosition);
  chain.SetBranchAddress("Track_Hit_ChargeSummedADC", &tTrack_Hit_ChargeSummedADC);
  chain.SetBranchAddress("Track_Hit_ChargeIntegral", &tTrack_Hit_ChargeIntegral);
  }
  if(version == "Feb"){
    chain.SetBranchAddress("NumberOfTracks_pmtrack",&tNumberOfTracks);
    chain.SetBranchAddress("Track_NumberOfHitsPerView_pmtrack",&tTrack_NumberOfHitsPerView);
    chain.SetBranchAddress("Track_Hit_Y_pmtrack", &tTrack_Hit_Y);
    chain.SetBranchAddress("Track_Hit_Z_pmtrack", &tTrack_Hit_Z);
  chain.SetBranchAddress("Track_Hit_dx_LocalTrackDirection_pmtrack", &tTrack_dx_LocalTrackDirection);
  chain.SetBranchAddress("Track_Hit_dx_3DPosition_pmtrack", &tTrack_dx_3DPosition);
  chain.SetBranchAddress("Track_Hit_ChargeSummedADC_pmtrack", &tTrack_Hit_ChargeSummedADC);
  chain.SetBranchAddress("Track_Hit_ChargeIntegral_pmtrack", &tTrack_Hit_ChargeIntegral);
  }
  for(int i=0; i<NEntries; i++){
    chain.GetEntry(i);
    int a=0;
    for(int j=0; j<tNumberOfTracks; j++){
      for(int k=0; k<2; k++){
        for(int l=a; l<a+tTrack_NumberOfHitsPerView[j][k]; l++){
//          int lem = find_lem(tTrack_Hit_Y[l], tTrack_Hit_Z[l]);
//          if(!isGood_lem(lem)){continue;}
//          if(!IsGoodChan(tTrack_Hit_Y[l], tTrack_Hit_Z[l])){continue;}
//          cout << tTrack_Hit_ChargeIntegral[l]/(ADC2CHARGE[k]*tTrack_dx_3DPosition[l]) << " " << tHit_GoodnessOfFit[l] << endl;
          
          HGoF[k].Fill(tTrack_Hit_ChargeIntegral[l]/(ADC2CHARGE[k]*tTrack_dx_3DPosition[l]),tHit_GoodnessOfFit[l]);
          
          HGoF_integral_local[k].Fill(tTrack_Hit_ChargeIntegral[l]/(ADC2CHARGE[k]*tTrack_dx_LocalTrackDirection[l]),tHit_GoodnessOfFit[l]);
          
//          HGoF_ByLEMs[find_lem(tTrack_Hit_Y[l], tTrack_Hit_Z[l])][k].Fill(tTrack_Hit_ChargeSummedADC[l]/(ADC2CHARGE[k]*tTrack_dx_3DPosition[l]),tHit_GoodnessOfFit[l]);
//          
//          HGoF_integral_local_ByLEMs[find_lem(tTrack_Hit_Y[l], tTrack_Hit_Z[l])][k].Fill(tTrack_Hit_ChargeIntegral[l]/(ADC2CHARGE[k]*tTrack_dx_LocalTrackDirection[l]),tHit_GoodnessOfFit[l]);
        }
        a+=tTrack_NumberOfHitsPerView[j][k];
      }
    }
  }
  string outpath = path_311data+"divers/";
  check_and_mkdir(outpath);
  outpath = outpath + to_string(run)+"/";
  check_and_mkdir(outpath);
  TFile ofile(string(outpath+"HGoF_vs_dqds.root").data(), "RECREATE");
  if(!ofile.IsOpen()){
    cout << " ERROR: TFile " << outpath << " can't be created " << endl;
    return;
  }
  #if verbose
  cout << " TFile " << outpath << "HGoF_vs_dqds.root has been created " << endl;
  #endif
  
  HGoF[0].Write();
  HGoF[1].Write();
  HGoF_integral_local[0].Write();
  HGoF_integral_local[1].Write();
  for( auto lem : lems ){
    HGoF_ByLEMs[lem][0].Write();
    HGoF_ByLEMs[lem][1].Write();
    HGoF_integral_local_ByLEMs[lem][0].Write();
    HGoF_integral_local_ByLEMs[lem][1].Write();
  }
  ofile.Close();
}//end macro
