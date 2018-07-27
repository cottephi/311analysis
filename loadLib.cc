#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

void loadLib(){

  //system("if [ -d ./obj ] ; then rm ./obj/* ; else mkidr ./obj ; fi");

  gStyle->SetPalette(kRainBow);
  gStyle->SetOptFit(11111);

  string currentdir = "/eos/user/p/pcotte/311analysis";
  std::string includepath="-I"+currentdir+"/lib";

  gSystem->SetBuildDir("/eos/user/p/pcotte/311analysis/obj",true);
  gSystem->AddIncludePath(includepath.c_str());
  gROOT->LoadMacro((currentdir+"/lib/311Lib.cc+").c_str());

  #define __INITIALIZED__
  
  //gApplication->Terminate();
  
  return;

}
