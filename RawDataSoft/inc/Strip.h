#ifndef STRIP_H_INCLUDED
#define STRIP_H_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>

using std::vector;
using std::string;
using std::cout;
using std::endl;

class Strip {
  public:
    // CONSTRUCTORS //
    Strip();
    Strip(int view, int strip_number, vector<int> adc_vec);
    ~Strip();

    // GETTERS //
    vector<int>  GetADCs() const { return m_adc_vec; }
    int GetTotalADC(){ return m_total_adc; }
    double GetMeanADC(){ return m_mean_adc; }
    double GetVarADC(){ return m_var_adc; }
    
    // SETTERS //
    void SetADCs(vector<int> adc_vec) { m_adc_vec = adc_vec; }
    void SetView(int view) { m_view = view; }
    void SetStripNumber(int strip_number) { m_strip_number = strip_number; }
    void AddADC(int adc) { m_adc_vec.push_back(adc); }

    // PRINTOUT //
    void PrintInfo() const;
    
    void Clear(){m_adc_vec.clear();}
    
    void ComputeTotalADC();

  private:
    int m_view;
    int m_strip_number;
    vector<int> m_adc_vec;
    int m_total_adc;
    double m_mean_adc;
    double m_var_adc;
};

#endif
