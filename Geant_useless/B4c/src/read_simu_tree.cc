#include <glob.h>
#include "read_simu_tree.hh"

using namespace std;

void Glob(const std::string& pattern, vector <string> &filenames) {
  // glob struct resides on the stack
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));

  // do the glob operation
  int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  if(return_value != 0) {
    globfree(&glob_result);
    stringstream ss;
    ss << "glob() failed with return_value " << return_value << endl;
    throw std::runtime_error(ss.str());
  }

  // collect all the filenames into a std::list<std::string>
  for(size_t i = 0; i < glob_result.gl_pathc; ++i) {
      filenames.push_back(string(glob_result.gl_pathv[i]));
  }

  // cleanup
  globfree(&glob_result);

  // done
  return;
}

void read_simu_tree(vector <MyPrimary> &primaries){

  string path = "~/Desktop/*.root";

  TChain chain("analysistree/anatree");
  vector <string> files;
  Glob(path,files);
  for(auto file : files){
    #if verbose
    cout << "Adding file: " << filename << endl;
    #endif
    chain.Add(filename.data());
  }
  #if verbose
  cout << "Files added. " << endl;
  #endif
  
  //read the tree and store all the tree information in the respective variables
  const int NMaxParticlesPerEvent=100;
  int NEventsPerRun=335;

  //Load ROOT file and tree
  int NEntries = (int)rTree->GetEntries();

  //Define variables to store the data of the ROOT file
  int MCTruth_Generator_NumberOfParticles;
  float MCTruth_Generator_StartTime[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartEnergy[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartMomentum_X[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartMomentum_Y[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartMomentum_Z[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartPoint_X[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartPoint_Y[NMaxParticlesPerEvent];
  float MCTruth_Generator_StartPoint_Z[NMaxParticlesPerEvent];
  

  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("MCTruth_Generator_NumberOfParticles",&MCTruth_Generator_NumberOfParticles);
  rTree->SetBranchAddress("MCTruth_Generator_StartTime",&MCTruth_Generator_StartTime);
  rTree->SetBranchAddress("MCTruth_Generator_StartEnergy",&MCTruth_Generator_StartEnergy);
  rTree->SetBranchAddress("MCTruth_Generator_StartMomentum_X",&MCTruth_Generator_StartMomentum_X);
  rTree->SetBranchAddress("MCTruth_Generator_StartMomentum_Y",&MCTruth_Generator_StartMomentum_Y);
  rTree->SetBranchAddress("MCTruth_Generator_StartMomentum_Z",&MCTruth_Generator_StartMomentum_Z);
  rTree->SetBranchAddress("MCTruth_Generator_StartPoint_X",&MCTruth_Generator_StartPoint_X);
  rTree->SetBranchAddress("MCTruth_Generator_StartPoint_Y",&MCTruth_Generator_StartPoint_Y);
  rTree->SetBranchAddress("MCTruth_Generator_StartPoint_Z",&MCTruth_Generator_StartPoint_Z);


  for(int i=0; i<NEntries; i++) //Event loop
  {
    rTree->GetEntry(i);

    //initialize classes (not pointers for the moment)

    for(int j=0; j<MCTruth_Generator_NumberOfParticles; j++) //Track loop
    {
      if(MCTruth_Generator_StartTime[j] == 0){
        MyPrimary dummy_MyPrimary;
        dummy_MyPrimary.start_E = MCTruth_Generator_StartEnergy[j];
        dummy_MyPrimary.start_x = MCTruth_Generator_StartPoint_X[j];
        dummy_MyPrimary.start_y = MCTruth_Generator_StartPoint_Y[j];
        dummy_MyPrimary.start_z = MCTruth_Generator_StartPoint_Z[j];
        dummy_MyPrimary.start_px = MCTruth_Generator_StartMomentum_X[j];
        dummy_MyPrimary.start_py = MCTruth_Generator_StartMomentum_Y[j];
        dummy_MyPrimary.start_pz = MCTruth_Generator_StartMomentum_Z[j];
        primaries.push_back(dummy_MyPrimary);
        continue;
      }
    }
  }
  return;
}
