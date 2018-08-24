#include "Event.h"

Event::Event(int event_number) : m_event_number(event_number)
{
//  m_total_adc = 0;
//  m_total_adcs[0] = 0;
//  m_total_adcs[1] = 0;
//  m_mean_adc = 0;
//  m_mean_adcs[0] = 0;
//  m_mean_adcs[1] = 0;
//  m_var_adc = 0;
//  m_var_adcs[0] = 0;
//  m_var_adcs[1] = 0;
  m_strips = {{},{}};
}

Event::Event(int event_number, vector<vector<Strip> > strips) : m_event_number(event_number), m_strips(strips)
{
//  m_total_adc = 0;
//  m_total_adcs[0] = 0;
//  m_total_adcs[1] = 0;
//  m_mean_adc = 0;
//  m_mean_adcs[0] = 0;
//  m_mean_adcs[1] = 0;
//  m_var_adc = 0;
//  m_var_adcs[0] = 0;
//  m_var_adcs[1] = 0;
}

Event::Event(int event_number, vector<Strip> strips0, vector<Strip> strips1) : m_event_number(event_number)
{
//  m_total_adc = 0;
//  m_total_adcs[0] = 0;
//  m_total_adcs[1] = 0;
//  m_mean_adc = 0;
//  m_mean_adcs[0] = 0;
//  m_mean_adcs[1] = 0;
//  m_var_adc = 0;
//  m_var_adcs[0] = 0;
//  m_var_adcs[1] = 0;
  m_strips = {strips0,strips1};
}

Event::~Event()
{

}

//void Event::ComputeTotalADC()
//{
//  for(size_t view = 0; view < m_strips.size(); view ++){
//    for(auto strip : m_strips[view]){
//      for(auto adc : strip.GetADCs()){
//        m_total_adcs[view] += adc;
//        m_total_adc += adc;
//      }
//    }
//  }
//  m_mean_adcs[0] = m_total_adcs[0]/320;
//  m_mean_adcs[1] = m_total_adcs[1]/960;
//  m_mean_adc = m_total_adc/1280;
//  for(size_t view = 0; view < m_strips.size(); view ++){
//    for(auto strip : m_strips[view]){
//      for(auto adc : strip.GetADCs()){
//        m_var_adcs[view] += pow((double)adc-m_mean_adcs[view],2);
//        m_var_adc += pow((double)adc-m_mean_adc,2);
//      }
//    }
//  }
//  m_var_ascs[0] = (double)m_var_ascs[0]/320;
//  m_var_ascs[1] = (double)m_var_ascs[1]/960;
//  m_var_asc = (double)m_var_asc/1280;
//}


void Event::Print() const
{
  cout << "    Event number: " << m_event_number << endl;
//  cout << "    Total ADC: " << m_total_adc << endl;
//  cout << "    Total ADC view 0: " << m_total_adcs[0] << endl;
//  cout << "    Mean ADC view 0: " << m_mean_adcs[0] << endl;
//  cout << "    Var ADC view 0: " << m_var_adcs[0] << endl;
//  cout << "    Total ADC view 1: " << m_total_adcs[1] << endl;
//  cout << "    Mean ADC view 1: " << m_mean_adcs[1] << endl;
//  cout << "    Var ADC view 1: " << m_var_adcs[1] << endl;
}
