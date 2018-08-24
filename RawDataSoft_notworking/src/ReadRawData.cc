#include "ReadRawData.h"


int main(int argc, char **argv){

//  TApplication theApp("App", &argc, argv);
  
  Strip strip;
  strip.SetStripNumber(2);
  Event event(840);
  cout << "chien" << endl;
  cout << strip.GetStripNumber() << " " << event.GetEventNumber() << endl;
  
//  TH2D h2("","",100,0,100,100,0,100);
//  h2.Fill(50,50,3);
//  h2.Draw("COLZ");

//  if(gPad) gPad->WaitPrimitive();
//  theApp.Terminate();
//  theApp.Run();
  
  return 0;
}
