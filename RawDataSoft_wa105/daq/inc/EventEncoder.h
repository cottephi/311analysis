////////////////////////////////////////////////////////////////////////
//
//
//  Class to produce raw data files / events  for testing purposes
//
//
////////////////////////////////////////////////////////////////////////

#ifndef __EVENTENCODER_H__
#define __EVENTENCODER_H__

#include <fstream>

#include "HuffDataCompressor.h"
#include "dlardaq.h"
 
namespace dlardaq
{
  class EventEncoder
  {
  public: 
    EventEncoder(size_t nch, size_t nsample);
    ~EventEncoder();
    
    // run number
    void SetRunNumber(uint32_t rn){ m_runnb = rn; }

    // compression flag
    void SetCompFlag(bool val){ m_dcflag = val; }

    // open input file
    int  Open(std::string finname, uint16_t flag);
    // filled with random data
    void WriteEvent( float ped, float rms );
    // filled with some actual data
    void WriteEvent(std::vector<adc16_t> &adc);
    void Close();
    
    // write to memory
    void WriteBufferEvent( float ped, float rms, std::vector<BYTE> &bytes );
    void WriteBufferEvent( std::vector<adc16_t> &adc,
			   std::vector<BYTE> &bytes );
    
    // in bytes
    size_t GetRunHeadSize() const { return m_RunHeadBytes; }
    size_t GetEventHeadSize() const { return m_EveHeadBytes; }
    size_t GetFooterSize() const { return m_EndFootBytes; }

    size_t GetTotEvents() const { return m_totev; }
    
  protected:
    void GenerateRandom(float ped, float rms, 
			std::vector<adc16_t> &adc);

    double genNoise(double mu, double sigma);

    void EncodeEvent( std::vector<adc16_t> &adc,
		      std::vector<BYTE> &buf );
    
    void WriteBytes(const BYTE* buf, size_t nb);
    
    BYTE* m_RunHeadBuf; // buffer to hold run header info
    BYTE* m_EveHeadBuf; // buffer to hold event header info
    BYTE* m_EndFootBuf; // buffer to hold file end info
    
    size_t m_RunHeadBytes;
    size_t m_EveHeadBytes;
    size_t m_EndFootBytes;

    void EncodeRunHead(uint16_t rflag);
    void EncodeEventHead(uint32_t evsize);
    void EncodeRunFoot();

    // total number of bytes written
    size_t m_totbytes;

    // basic data parameters
    short  m_nadc;
    size_t m_nch;
    size_t m_nsample;
    size_t m_datasz;

    // total number of events in the file
    uint32_t m_totev;
    uint32_t m_runnb;
    
    // data compression flag
    bool m_dcflag;
    
    // input file
    std::ofstream fOut;
  };

};


#endif
