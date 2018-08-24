#include "Strip.h"
Strip::Strip()
{
  m_view = -1;
  m_strip_number = -1;
  m_adc_vec = {};
  m_total_adc = 0;
  m_mean_adc = 0;
  m_var_adc = 0;
}
Strip::Strip(int view, int strip_number, vector<int> adc_vec) : m_view(view), m_strip_number(strip_number), m_adc_vec(adc_vec)
{
  m_total_adc = 0;
  m_mean_adc = 0;
  m_var_adc = 0;
}

Strip::~Strip()
{
  
}

void Strip::Print() const
{
  cout << "    Total ADC: " << m_total_adc << endl;
  cout << "    Mean ADC: " << m_mean_adc << endl;
  cout << "    Var ADC: " << m_var_adc << endl;
}

void Strip::ComputeTotalADC()
{
  for(auto adc : m_adc_vec){
    m_total_adc += adc;
  }
  m_mean_adc = (double)m_total_adc/1667;
  for(auto adc : m_adc_vec){
    m_var_adc += pow(((double)adc-m_mean_adc),2);
  }
  m_var_adc = (double)m_var_adc/1667;
}
