////////////////////////////////////////////////////////////////////////////////
//
// 311 common analysis function and variables
//
////////////////////////////////////////////////////////////////////////////////


#ifndef  __311LIB_H
#define __311LIB_H

//#ifdef __MAKECLING__
#pragma link C++ class track+;
//#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TMultiGraph.h>

#define verbose 1


using namespace std;

//Default paths and values************************************************************

string path_311data;
string path_wa105_311data;
double lem_size;
vector<int> tpc_boundaries;
double pitch;
double dQdx_cut_min;
double dx,dy,dz;
bool IsBatch = false;

//cuts********************************************************************************
double length_cut = 0;
double theta_cut = -1.;
double phi_cut = -1.;
double ds_cut = 1e9;
double GoF_cut = 1e9;
bool only_throughgoing = false;
bool only_throughgoing_x = false;
bool only_throughgoing_y = false;
bool only_throughgoing_z = false;

double rho_ref = 11.2644; //980/87, from 3L
double rho_run;
double rho_var_run;
TH2D h_ExtrEff_vs_LemExtr;
TH2D h_IndEff_vs_LemInd;
TGraph gr_rho;
TGraphErrors Gain_graph_3L;
TH1D h_rho;
vector<double> simulated_extr;
vector<double> simulated_amp;
vector<double> simulated_ind;
map<int,double> gain_corrections;
vector<int> tracks_selected_by_highway;
map<int, map<string,float> > runs_and_fields;
pair<double,double> Arho_from_3L = make_pair(-1,-1); //pressure 980mbar, temperature 87K
pair<double,double> Brho_from_3L = make_pair(-1,-1);

const string path_MC_input = "/eos/experiment/wa105/offline/LArSoft/MC/MC6/ROOT/g4detsim/";
const string path_311analysis = "/eos/user/p/pcotte/311analysis/";
const string path_MC_output = "/eos/user/p/pcotte/311data/MC/";
const string dqds_cosmics = path_311analysis + "dqds_cosmics.dat";
const string runs_file = path_311analysis + "all_runs.csv";
const string runs_headers = "/eos/user/p/pcotte/311data/runs_headers/";
const string slow_control = "/eos/user/p/pcotte/311data/slow_control/";
const string highway_output = "/eos/user/p/pcotte/from_github_311analysis/Analysis/Event-track-selection/HighwayAlgorithm/HighwayOutput/text_files/";

string SelectTrack_Input;
string SelectTrack_Output;
string SelectTrack_MC_Output;
string YZ_SelectTrack_Output;
string GainAna_Output;
string GainAna_dQds_Output;
string PlotGainAna_Input;
string PlotGainAna_Output;
string Purity_Output;
string Dodqdx_Output;
string cuts_analysis_Output;
string MPV_vs_stuff_Output;
string MPV_vs_stuff_tracks_Output;
string gain_by_lem_Output;
string YZ_Output;
string dQds_Output;
string dQds_charging_up_Output;
string dQds_YZ_Output;
string dQds_MC_Output;
string charging_up_Output;
string gain_stability_Output;

vector<int> MyColorPalette = {632+2,632,800+9,800+7,800-3,800-2,820+10,820-9,820,840+8,840+7,840,860+8,860+7,860,880+7,880+1,880};
vector<string> params = {"date","TE0037","TE0038","TE0039","TE0040","TE0041","TE0042","TE0043","TE0044","TE0045","TE0046","TE0047","TE0048","TE1001","TE1002","TE1003","TE1004","TE0049","TE0050","TE0051","TE0052","TE0053","TE0054","TE0055","TE0056","TE0057","TE0058","TE0059","TE0060","PE0006"};

float ExtrField_Gushin[23] = {0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
float ExtrEff_Gushin[23] = {0.3, 0.35, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.85, 0.9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
TGraph GushinEff(23, ExtrField_Gushin, ExtrEff_Gushin);
TH1F h_GushinEff("h_GushinEff","h_GushinEff",23,0.7,3.0);

//Geometry**********************************************************************

bool Load_Version(string version);
bool load_runs_2D(string scan_type, int scan_num, vector<int> &run_list, vector<pair<float, float> > &field, vector<string> &const_fields_names, vector<float> &const_fields_values);
bool load_runs(string scan_type, int scan_num, vector<int> &run_list, vector<float> &field, vector<string> &const_fields_names, vector<float> &const_fields_values);
bool load_cosmics();
bool load_run_lists();
bool load_eff_simu_graphs();
bool load_rho_run(int run);
void load_fit_3L();
void load_Gushin_Eff();
void load_gain_eff_corrections();
bool load_highway(int run);

const int NUM_OF_VIEWS = 2;
const int NUM_OF_LEMS = 12;
const int NUM_OF_GOOD_LEMS = 8;
double N_Ch_0 = 320;
double N_Ch_1 = 960;
int tdc = 1667;
/*double ADC2CHARGE = 45.31875; //ADC*ticks (from qScan)*/
double ADC2CHARGE[2] = {55, 67};
vector<int> bad_channels = {576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607};
vector<int> lems = {2, 4, 5, 6, 7, 8, 9, 11}; //active lems
vector<int> vol_cut = {0, 0, 0, 0, 0, 0}; //in cm
static int default_int = 1;
double mpv_cosmics = -1.;
const double drift_velocity = 150. ; //cm/ms
double T0 = 87;
double P0=980;
double corr = 1.;
double e_lifetime_const = 7.;//ms //Caspar talk 29 June
string method_ds = "3D";
string method_dQ = "sum";
bool highway = false;

double golden_ratio = 1.618;

//Data products ****************************************************************

class hit{
  public:
    int run;
    int subrun;
    int event;
    int view;
    int track_id;
    double peak_amp;
    double peak_time;
    double chi2;
    double sp_x;
    double sp_y;
    double sp_z;
    double dq_sum;
    double dq_integral;
    double dx_3D;
    double dx_local;
    double dx_MC;
    double purity_correction;
    int lem;
    double gain_density_correction_factor;
    double GoF;
};

class track{
  public:
    int run;
    int subrun;
    int event;
    int id;
    double start_x;
    double start_y;
    double start_z;
    double start_theta;
    double start_phi;
    double start_vx;
    double start_vy;
    double start_vz;
    double end_x;
    double end_y;
    double end_z;
    double end_theta;
    double end_phi;
    double end_vx;
    double end_vy;
    double end_vz;
    double length;
    double theta;
    double phi;
    vector<double> MPVs_dQds;
    int nhits;
//    vector<free_hit> free_hits_trk;
    vector<hit> hits_trk;
};

//#ifdef __MAKECLING__ 
#pragma link C++ class vector<track>+; 
//#endif

//Functions declaration*********************************************************

int find_lem(double y, double z);

pair<int,int> find_channels(double y, double z);

vector<int> find_YZ(int lem);

bool isGood_lem( int lem );

bool IsGoodChan(int ch0, int ch1);
bool IsGoodChan(double y, double z);
bool IsGoodChan(hit myhit);

double find_projection(hit h);

void read_tree_Feb(TChain *rTree, vector<track> & tracks, int &tstart, int &tend, int to_read = 0);
void read_tree_June(TChain *rTree, vector<track> & tracks, int &tstart, int &tend, int to_read = 0);

int drays_mitigation(track & t);

double get_theta(track t);

double get_phi(track t);

vector<double> get_dss(track t);

vector<double> get_dss_from_angles(double theta, double phi);

double find_closest_fit_value_in_bin(TH1D* myhisto, int bin, TF1 *myfunc);

double homemade_get_pvalue(TF1 *myfunc, TH1D *myhisto);
double truncated_mean(std::vector<double> &vec, double rmlow, double rmhigh );

void FWHM(double &st, double &end, TH1D *h);

double find_max_bin(TH1D *hist, int min_bin = 3);

double langaufun(double *x, double *par);

TF1 *langaufit(TH1D *his, double *fitrange, vector<double> startvalues,
 vector<double> parlimitslo, vector<double> parlimitshi, vector<double> &fitparams, vector<double> &fiterrors, double &pvalue, double pvaluelim, int &NDF, bool find_best = false, bool gauss = false);

void plot_selected_tracks(vector<track> tracks, string path, int ntracks = -1);

////////////not 311-related/////////////

inline vector<string> glob(const string& pat);

inline bool ExistTest (const std::string& name);

vector<double> fit_dQds(TH1D *hdQds, bool gauss = false, int min_number_of_hits = 500, double pvaluelim = 0.05, double sigmalim = 1000, TGraphErrors* graph = 0, string scan_type = "", float x = 0, float xer = 0);

void fill_2d_graph(TGraph2D* graph, double x, double y, double z);

void fill_2d_hist(TH2D* h, double x, double y, double z);

vector<double> ReadFit(TH1D* hdQds, TGraphErrors* graph = 0, string scan_type = "", float x = 0, float xer = 0);

bool init_graph_gain(vector<TGraphErrors*> &mpv_field, map<int, vector<TGraphErrors*> > &mpv_field_ByLEMs, vector<TMultiGraph*> &mpv_field_AllLEMs, vector<int> &scan_nums_for_AllLEMs, string scan_type = "", int scan_num = 0);

bool init_graph_purity_and_charging_up(vector<TGraphErrors> &graph, map<int, vector<TGraphErrors>> &graph_ByLEMs, string name = "");

void sum_views(vector<TGraphErrors> &graph, map<int, vector<TGraphErrors> > &graph_ByLEMs);
void sum_views(vector<TGraphErrors> &graph, map<int, vector<TGraphErrors*> > &graph_ByLEMs);

void sum_views_multigraph(vector<TMultiGraph*> &graph);

void sum_views_singlegraph(vector<TGraphErrors> &graph);
void sum_views_singlegraph(vector<TGraphErrors*> &graph);

double e_lifetime(TGraphErrors *graph);

double RFromBirk(double drift, double MPV = -1.);

double MPVs_from_theory(double ds, double drift = 0.5);

double get_eff(double extraction = -1., double amplification = -1., double induction = -1.);

pair<TH2D*, TGraph*> normalise_gain_histo_2d(TH2D* gr, string scan_type);

void get_eff_multigraphs(TMultiGraph *mg, string scan_type, vector<int> scan_nums, TFile *ofile, string eff_to_use = "simu", string draw = "");

void get_eff_graphs(TGraphErrors *gr, string scan_type, int scan_num, string eff_to_use = "simu", string draw = "");

TMultiGraph* normalise_gain_graph(TMultiGraph *mg, string scan_type, vector<string> effs_to_use);

double GetMaxOFMapMapTH1D(map<double,map<int,TH1D> > MapMapTH1D);

void draw_gain_dQds(map<double,map<int,TH1D> > mpv_field, string scan_type, int scan_num, TFile *ofile, string path = "", int ilem = -1);

void draw_gain_multigraph(TMultiGraph *mg, string scan_type, TFile *ofile, bool draw_leg = false, string draw = "");

void draw_gain_graph(TGraphErrors *gr, string scan_type, int scan_num, TFile *ofile, string draw = "");

void draw_hist_2d(TH2D* h, string scan_type, TFile *ofile, string directory);

void fit_gain(TGraphErrors* gr, vector<double> par);

double correct_dx_for_lifetime(double dx, double e_lifetime = e_lifetime_const);

void fit_charging_up(TGraphErrors &gr);

double theoretical_gain(double T, double rho, double E);

double gain_correction_for_rho(double rho, double E);

double get_density_for_hit(int hit_time, int run);

int read_or_do_fit(string filename_fitted, string filename_nonfitted, bool recreate_fit_file, string &ifilename, int &files_not_found);

bool init_histo_track_cuts(vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected, int time = -1);

void rec_track_dQds(track t, vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected, map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected);

bool select_tracks(string cut_type, vector<track> tracks, vector<track> & mips, vector<TH1D> &hdQds, map<int, vector<TH1D> > &hdQds_ByLEMs, map<int, vector<TH1D> > &hdQds_ByDx, map<int, map<int, vector<TH1D> > > &hdQds_ByDx_ByLEMs, vector<TH1D> &hdQds_Dx_Corrected,  map<int, vector<TH1D> > &hdQds_ByLEMs_Dx_Corrected);

void check_and_mkdir(string path);

bool get_histo_in_inputfile(TH1D &hdQds, TFile *runfile, string name_to_get, bool &read_fit);

string set_cuts(string cut_type);

bool pressure(string srun = "801");

#endif // __311LIB_H
