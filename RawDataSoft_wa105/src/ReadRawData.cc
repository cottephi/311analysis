#include "ReadRawData.h"



void PrintHelp(){
  cout << "options to use: " << endl;
  cout << "  -r run_number (requiered)" << endl;
  cout << "  -s subrun_number. If not specified, will read all the subruns" << endl;
  cout << "  -f first event to read" << endl;
  cout << "  -l last event to read. If identical to first event, only one will be read. If inferior to first event, they are switched." << endl;
  cout << "If only first is specified, will read from it to the end. If only last is specified, will read from start to it." << endl;
  cout << "IMPORTANT: event numbers are understood as global. So subruns have n events, the first event of the second subrun is n+1." << endl;
  cout << "If the last event number is superior to total event in run, the loop will just end at the last event without errors." << endl;
  return;
}

int main(int argc, char **argv){
  int global_counter = 0;
  int run = -1; 
  int subrun = -1;
  int first_event = -1;
  int last_event = -1;
  char cc;

  while((cc = getopt(argc, argv, "r:s:f:l:h")) != -1){  
    switch(cc){
      case 'r': 
        run =  atoi(optarg);
        break;  

      case 's':
        subrun = atoi(optarg); 
        break;
         
      case 'f':
        first_event = atoi(optarg);  
        break; 

      case 'l':
        last_event = atoi(optarg);
        break; 

      case 'h':
        PrintHelp();
        exit(1);  
    }
  }
  if(argc < 2 or run < 0){
    PrintHelp();
    return 1;
  }
  string rawdata_folder = "/eos/experiment/wa105/data/311/rawdata/" + to_string(run) + "/";
  if(!ExistTest(rawdata_folder)){
    cout << "Directory " << rawdata_folder << " not found." << endl;
    return 1;
  }
   
  vector<string> input_files;

  if(subrun < 0 and first_event < 0 and last_event < 0){
    first_event = 0;
    last_event = 1e9;
    cout << "Processing all subruns and events in " << rawdata_folder << "..." << endl;
    string wildcard_path = rawdata_folder + "*";
    for( auto ssubrun : glob(wildcard_path) ){
      cout << "  Adding file " << ssubrun << endl;
      input_files.push_back(ssubrun);
    }
  }
  else if(subrun > 0 and first_event < 0 and last_event < 0){
    first_event = 0;
    last_event = 1e9;
    string file = rawdata_folder + to_string(run) + "-" + to_string(subrun) + ".dat";
    if(!ExistTest(file)){
      cout << "File " << file << " not found" << endl;
      return 1;
    }
    cout << "Processing all events in subrun " << subrun << " of run " << run << endl;
    input_files.push_back(file);
  }
  else if(subrun < 0 and (first_event >= 0 or last_event >= 0) ){
    if(first_event < 0){
      first_event = 0;
      subrun = SubRunOfEvent(last_event);
      string wildcard_path = rawdata_folder + "*";
      int isubrun = -1;
      cout << "Processing from first event to event " << last_event << " (subrun) " << subrun << endl;
      for( auto ssubrun : glob(wildcard_path) ){
        isubrun = atoi(ssubrun.substr(ssubrun.find_first_of("-")+1,ssubrun.find_first_of(".")-ssubrun.find_first_of("-")-1).data());
        if(isubrun > subrun){continue;}
        cout << "  Adding subrun " << isubrun << endl;
        input_files.push_back(ssubrun);
      }
    }
    else if(last_event < 0){
      last_event = 1e9;
      subrun = SubRunOfEvent(first_event);
      string wildcard_path = rawdata_folder + "*";
      int isubrun = -1;
      cout << "Processing from event " << first_event << " (subrun) " << subrun << " to last event" << endl;
      for( auto ssubrun : glob(wildcard_path) ){
        isubrun = atoi(ssubrun.substr(ssubrun.find_first_of("-")+1,ssubrun.find_first_of(".")-ssubrun.find_first_of("-")-1).data());
        if(isubrun < subrun){continue;}
        cout << "  Adding subrun " << isubrun << endl;
        input_files.push_back(ssubrun);
      }
    }
    else if(first_event >= 0 and last_event >= 0){
      if(first_event > last_event){
        int dummy = first_event;
        first_event = last_event;
        last_event = dummy;
      }
      int first_subrun = SubRunOfEvent(first_event);
      int last_subrun = SubRunOfEvent(last_event);
      string wildcard_path = rawdata_folder + "*";
      int isubrun = -1;
      cout << "Processing from event " << first_event << " (subrun) " << first_subrun << " to event " << last_event << " (subrun) " << last_subrun << endl;
      for( auto ssubrun : glob(wildcard_path) ){
        isubrun = atoi(ssubrun.substr(ssubrun.find_first_of("-")+1,ssubrun.find_first_of(".")-ssubrun.find_first_of("-")-1).data());
        if(isubrun < first_subrun or isubrun > last_subrun){continue;}
        cout << "  Adding subrun " << isubrun << endl;
        input_files.push_back(ssubrun);
      }
    }
  }
  else{
    cout << "Cannot specify event numbers and subrun at the same time!" << endl;
    return 1;
  }
  if(input_files.size() == 0){
    cout << "No subruns found" << endl;
    return 1;
  }
  
  
//  vector<Event> events = {};
  for(auto ssubrun : input_files){
    int isubrun = atoi(ssubrun.substr(ssubrun.find_first_of("-")+1,ssubrun.find_first_of(".")-ssubrun.find_first_of("-")-1).data());
    cout << "  Reading subrun " << isubrun << endl;
    dlardaq::EventDecoder DaqDecoder( 1280, 1667 ); 
    DaqDecoder.Open(rawdata_folder+ssubrun);
    size_t NEvent = DaqDecoder.GetTotEvents();
    for(size_t event = 0; event < NEvent; event++){
    	global_counter++;
    	if(global_counter > 10){return 0;}
      int glob_event = event + 335*isubrun;
      if(glob_event < first_event or glob_event > last_event){continue;}
  		string path_outfile = "/eos/user/p/pcotte/311analysis/RawDataSoft/textfiles/" + to_string(run) + "/" + to_string(glob_event) + ".txt";
	    cout << "    Saving event " << glob_event << " in " << path_outfile << endl;
  		check_and_mkdir(path_outfile);
  		std::ofstream outfile(path_outfile, std::ofstream::out);
      dlardaq::evheader_t event_header; 
      vector<dlardaq::adc16_t> adc_vector;
      DaqDecoder.GetEvent( event, event_header, adc_vector);
//      Event myevent(event);
//      Strip strip;
      for(size_t i = 0; i < adc_vector.size(); i++){
        if(i%1667 == 0 and i > 0){
          outfile << endl;
//          int strip_number = i/1667;
//            strip.SetStripNumber(i/1667);
//          if(strip_number < 320){
//            strip.SetView(0);
//            myevent.AddStrip(0,strip);
//          }
//          else if(strip_number < 1280){
//            strip.SetView(1);
//            myevent.AddStrip(1,strip);
//          }
//          else{
//            cout << "ERROR: strip numbering inconsistent (strip number is " << strip_number << ")" << endl;
//            return 1;
//          }
//          strip.Clear();
        }
        outfile << adc_vector[i] << " ";
//        strip.AddADC(adc_vector[i]);
      }
//      myevent.PrintInfo();
//      cout << "    Event size: " << event_header.ev_size/(1024*1024) << endl;
//      cout << "    Event quality flag: " << bitset<8>(event_header.dq_flag) << endl;
//      cout << "    Trigger number: " << event_header.trig_info.num << endl;
//      cout << "    Event timestamp: " << event_header.trig_info.ts.tv_sec <<" s " << event_header.trig_info.ts.tv_nsec << " ns" << endl;
//      events.push_back(myevent);
      outfile.close();
    }
  }
//  if(events.size() == 0){
//    cout << "No events analysed" << endl;
//    return 1;
//  }
  
//  int i = 0;
//  for(auto event : events){
//    vector<TH2D> h_event_display = {};
//    h_event_display.push_back(TH2D("",string("View 0, event " + to_string(event.GetEventNumber())).data(),320,0,320,1667,0,1667));
//    h_event_display.push_back(TH2D("",string("View 1, event " + to_string(event.GetEventNumber())).data(),960,0,960,1667,0,1667));
//    for(auto strips : event.GetStrips()){
//      for(int strip = 0; strip < strips.size(); strip++){
//        for(int iadc = 0; iadc < strips[strip].GetADCs().size(); iadc++){
//          h_event_display[strips[strip].GetView()].Fill(strip,iadc,strips[strip].GetADCs()[iadc]);
//        }
//      }
//    }
//    i++;
//    if(i == 10){
//      
//      break;
//    }
//  }
  
  
  return 0;
}
