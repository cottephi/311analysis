#ifndef READRAWDATA_H_INCLUDED
#define READRAWDATA_H_INCLUDED

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <glob.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <unistd.h>
#include <bitset>
#include <sstream>
#include <iterator>

//#include <TCanvas.h>
//#include "Event.h"
//#include "TROOT.h"
//#include <TH2D.h>
#include "EventDecoder.h"


using std::vector;
using std::string;
using std::bitset;
using std::cout;
using std::endl;
using std::istringstream;
using std::istream_iterator;
using std::to_string;

inline bool ExistTest (const string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

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

int SubRunOfEvent(int event){return event/335;}

class WordDelimitedBySlash : public std::string
{};
std::istream& operator>>(std::istream& is, WordDelimitedBySlash& output)
{
   std::getline(is, output, '/');
   return is;
}

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
      cout << "Creating directory " << partpath << endl;
      mkdir(partpath.data(),0777);
    }
    partpath += "/";
  }
  return true;
}



#endif
