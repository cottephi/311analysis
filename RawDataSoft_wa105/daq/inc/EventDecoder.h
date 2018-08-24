//////////////////////////////////////////////////////////////////////////////////////
//
//  Decoder for the DAQ events for charg readout
//  The general idea is to try to make this object threadsave when possible
//
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef __EVENTDECODER_H__
#define __EVENTDECODER_H__

// for mutexes
//#include <pthread.h>

#include "dlardaq.h" // in daq folder
#include <fstream>

// for compressed data
#include "HuffDataCompressor.h" // in daq folder

// max time to wait to receive data 
#define TMAXWAIT 2

namespace dlardaq
{
  class EventDecoder
  {
  public: 
    EventDecoder();
    EventDecoder(size_t nch, size_t nsample);
    ~EventDecoder();

    //void SetNCh(size_t val){ m_nch = val; }
    //void SetNSample(size_t val){ m_nsample = val; }

    size_t GetNCh() const { return m_di.nch; }
    size_t GetNSample() const { return m_di.nsample; }
    
    // open input file
    ssize_t Open(std::string finname);
    void Close();
    
    // get a given event from fil
    ssize_t GetEvent( size_t evnum, dlardaq::evheader_t &eh, 
		      std::vector<adc16_t> &adc );

    // get event from buffer
    ssize_t GetEvent( dlardaq::evheader_t &eh, 
		      std::vector<adc16_t> &adc );
    
    // read event from online buffer
    // and store data internally
    void ReadBuffer(const char *buf, size_t nb);

    
    //
    bool Compressed() const { return m_dcflag; }

    
    size_t GetTotEvents() const { return m_totev; }

    runheader_t GetRunHeader()   { return m_rnh; }// runheader_t defined in dlardaq.h
    evheader_t  GetEventHeader() { return m_evh; }// evheader_t defined in dlardaq.h
    footer_t    GetFileFooter()  { return m_flf; }// footer_t defined in dlardaq.h
    datainfo_t  GetDataInfo() { return m_di; }// datainfo_t defined in dlardaq.h
    


  protected:
    void Init();

//    int lock(pthread_mutex_t &mutex) { return  pthread_mutex_lock(&mutex); }
//    int trylock(pthread_mutex_t &mutex) { return  pthread_mutex_trylock(&mutex); }
//    int unlock(pthread_mutex_t &mutex) { return  pthread_mutex_unlock(&mutex); }   
    

    // decode event bytes
    void ReadEvent(std::vector<adc16_t> &adc, bool headonly = false);
    
    // read a byte fector from current position
    void ReadBytes( std::vector<BYTE> &bytes );
    ssize_t Decode( const char *buf, size_t nb, bool cflag,
		    std::vector<adc16_t> &adc);
    
    bool IsFirstPacket(const char *buf, size_t nb);

    //
    int m_Verbosity;  //

    std::vector<BYTE> m_RunHeadBuf; // BYTE definded in dlardaq.h. buffer to hold run header info
    std::vector<BYTE> m_EveHeadBuf; // buffer to hold event header info
    std::vector<BYTE> m_EndFootBuf; // buffer to hold file end info
    std::vector<BYTE> m_EveDataBuf; // event data buffer

    std::vector<adc16_t> m_EveData; // adc16_t definded in dlardaq.h. event adc data
    
    //
    ssize_t m_bytes_left;
    
    runheader_t m_rnh;// runheader_t defined in dlardaq.h
    evheader_t  m_evh;// evheader_t defined in dlardaq.h
    footer_t    m_flf;// footer_t defined in dlardaq.h

    //
    std::streampos m_pstart;
    std::streampos m_pend;
    std::vector<std::streampos> m_events;
    
    // total number of events in the file
    size_t m_totev;
    
    // basic data parameters
    //short  m_nadc;
    //size_t m_nch;
    //size_t m_nsample;
    datainfo_t  m_di; //datainfo_t defined in dlardaq.h

    // data compression flag
    bool m_dcflag;
    
    // input file
    std::ifstream m_file;

    // mutex
//    pthread_mutex_t m_data_mutex;
  };
};

#endif
