/////////////////////////////////////////////////////////////////////////////////
//
//  declarations / primitives / functions to handle raw data from DAQ
//  ADC resolution is 12bit
// 
//  Created: vgalymov, Sat Jul  2 15:25:46 CEST 2016
//  Modified: 
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __DLARDAQ_H__
#define __DLARDAQ_H__

#include <exception>
#include <vector>
#include <string>
#include <ctime>
#include <stdint.h>

// key for raw data
#define EVSKEY 0xFF
#define ENDKEY 0xF0

// bit check
#define SETBYTEBIT(var, pos) ( var |= 1 << pos )
#define CLEARBYTEBIT(var, pos) ( var &= ~(1 << pos) )
#define CHECKBYTEBIT(var, pos) ( (var) & (1<<pos) )

// flag bit for compression
#define DCBITFLAG 0x6 // 0x0 LSB -> 0x7 MSB
#define GETDCFLAG(info) (CHECKBYTEBIT(info, DCBITFLAG)>0)
#define SETDCFLAG(info) (SETBYTEBIT(info, DCBITFLAG))

// event data quality flag
#define EVDQFLAG(info) ( (info & 0x3F ) == 0 )

///
namespace dlardaq
{
  typedef char BYTE;
  typedef uint16_t adc16_t;

  //
  // formatting exception
  class formatexception : public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Bad file format";
    }
  };
  
  static formatexception fex;// in dlardaq.h
  
  // specify sizes of different headers
  // run header in bytes
  static const size_t RunHeadSz  = 5;
  // event header in bytes
  static const size_t EveHeadSz  = 35;
  // run footer in bytes
  static const size_t FileFootSz = 4;
  
  // adc bits
  static const short  BitsADC = 12;
  
  // trigger structure from WR trig handler
  typedef struct trigger_t
  {
    uint8_t type;
    uint32_t num;
    struct timespec ts; //{ time_t ts.tv_sec, long tv_nsec }
  } trigger_t;
 
  // structure to hold decoded run header
  typedef struct runheader_t
  {
    uint32_t run_num;
    uint8_t  run_flags; 
  } runheader_t;

  // detector data information
  typedef struct datainfo_t
  {
    std::string name;
    short bitadc;    // ADC resolution
    size_t nch;      // total number of channels
    size_t nsample;  // total ADC samples per ch
  } datainfo_t;


  // strucutre to hold decoded event header
  typedef struct evheader_t
  {
    trigger_t trig_info; // trigger info
    uint8_t  dq_flag;    // data quality flag
    uint32_t ev_num;     // event number
    uint32_t ev_size;    // size of event in bytes
  } evheader_t;

  // structure to hold footer information
  typedef struct footer_t
  {
    uint16_t num_events;
  } footer_t;
  
  // pack 16 bits 12 bits 
  void pack16into12(const void *in, void *out, size_t n);
  
  // unpack 12 bits into 16 bits
  void unpack12into16(const void *in, void *out, size_t n);
  // pack 16 bits 12 bit

  // write binary file with 12 bit words
  void write12(const char *fname, std::vector<adc16_t> &adc);
  // read binary file from 12 bit  words
  void read12(const char *fname, std::vector<adc16_t> &adc);
  

  // get byte content for a given data type
  // NOTE: cast assumes host byte order
  template<typename T> T ConvertToValue(const void *in)
    {
      const T *ptr = static_cast<const T*>(in);
      return *ptr;
    }

  // functions to decode header informations from binary stream
  ssize_t decode_runhead(const char* buf, runheader_t &rh);
  ssize_t decode_evehead(const char* buf, evheader_t &eh);
  ssize_t decode_filefoot(const char* buf, footer_t &rf);

  // get basic detector information from the decoded run header
  ssize_t get_datainfo( const runheader_t &rh, datainfo_t &di);


  // properties of raw data produced by DAQ
  namespace raw311
    {
      static const std::string name = "3x1x1";
      static const size_t nsa       = 1667; // total number of samples per ch
      static const short nchcard    = 64;   // number of channels per card
      static const short nccrate    = 5;    // cards per crate
      static const short ncrate     = 4;    // crates

      // total number of channels
      static const size_t nch = 1280; // == nchcard *nccrate * ncrate; 
    }

  namespace raw666
    {
      static const std::string name = "6x6x6";
      static const size_t nsa       = 10000; // total number of samples per ch
      static const short nchcard    = 64;    // number of channels per card
      static const short nccrate    = 10;    // cards per crate
      static const short ncrate     = 12;    // crates

      // total number of channels
      static const size_t nch = 7680; // == nchcard *nccrate * ncrate; 
    }
};

#endif
