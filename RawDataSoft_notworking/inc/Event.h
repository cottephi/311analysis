#ifndef EVENT_H_INCLUDED
#define EVENT_H_INCLUDED


#include <vector>
#include <iostream>

#include "Strip.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;

class Event{
  public:
    // CONSTRUCTORS //
    Event(int event_number);
    Event(int event_number, vector<vector<Strip> > strips);
    Event(int event_number, vector<Strip> strips0, vector<Strip> strips1);
    ~Event();

    // GETTERS //
    int  GetEventNumber() const { return m_event_number; }
    vector<vector<Strip> > GetStrips() const { return m_strips; }
    vector<Strip> GetStrips(int view) const { return m_strips[view]; }
//    int GetTotalADC() { return m_total_adc; }
    
    // SETTERS //
    void AddStrip(int view, Strip strip) { m_strips[view].push_back(strip); }
    void SetEventNumber(int event_number) { m_event_number = event_number; }
    void SetStrips(int event_number, vector<vector<Strip> > strips) { m_strips = strips; }
    void SetStrips(int event_number, vector<Strip> strip0, vector<Strip> strip1) { m_strips = {strip0, strip1}; }
    
//    void ComputeTotalADC();

    // PRINTOUT //
    void Print() const;
    
  private:
    float  m_event_number; // 
    vector<vector<Strip> > m_strips;
//    int m_total_adc;
//    int m_total_adcs[2];
//    int m_mean_adc;
//    int m_mean_adcs[2];
//    int m_var_adc;
//    int m_var_adcs[2];

};
#pragma link C++ class vector<Event>+; 

#endif
