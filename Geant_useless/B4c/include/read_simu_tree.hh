#ifndef  __READ_SIMU_TREE_H
#define __READ_SIMU_TREE_H

#pragma link C++ class MyPrimary+;

#define verbose 1

class MyPrimary{
  public:
    double start_E;
    double start_x;
    double start_y;
    double start_z;
    double start_px;
    double start_py;
    double start_pz;
};

class MyHit{
  public:
    int view;
    double dE;
    double x;
    double y;
    double z;
    double ds;
};

#pragma link C++ class vector<MyPrimary>+; 
