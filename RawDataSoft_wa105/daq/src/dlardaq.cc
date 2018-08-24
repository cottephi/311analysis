/////////////////////////////////////////////////////////////////////////////////
//
//  declarations / primitives / functions to handle raw data from DAQ
//  ADC resolution is 12bit
//
//  Created: vgalymov, Sat Jul  2 15:25:46 CEST 2016
//  Modifed: 
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <string.h>
#include <assert.h>

#include "dlardaq.h"
//#include "LogMsg.h"

using namespace std;
using namespace dlardaq;

//
// n8bit is the size of the output that should be 
// allocated to be x3/2 larger of the input
//
void dlardaq::pack16into12(const void *in, void *out, size_t n8bit)
{
  const adc16_t* in16 = static_cast<const adc16_t*>(in);
  BYTE* out8 = static_cast<BYTE*>(out);// in dlardaq.h
  BYTE* stop = out8 + n8bit;// in dlardaq.h
  
  while(out8!=stop)
    {
      // we pack 2x12 bit values into 3x8bit
      adc16_t v1 = *in16++;// in dlardaq.h
      adc16_t v2 = *in16++;// in dlardaq.h
      
      *out8++ = ( v1 >> 4 ) & 0xff;
      *out8++ = ( (v1 & 0xf) << 4 ) + ( ( v2 >> 8) & 0xf );
      *out8++ = v2 & 0xff;
    }
}

//
// n8bit is the size ouf the input that should be x3/2
// larger than the output vector
//
void dlardaq::unpack12into16(const void *in, void *out, size_t n8bit)
{
  const BYTE* in8  = static_cast<const BYTE*>(in);// in dlardaq.h
  const BYTE* stop = in8 + n8bit;// in dlardaq.h

  adc16_t* out16 = static_cast<adc16_t*>(out);// in dlardaq.h

  while(in8!=stop)
    {
      BYTE v1 = *in8++;// in dlardaq.h
      BYTE v2 = *in8++;// in dlardaq.h
      BYTE v3 = *in8++;// in dlardaq.h

      *out16++ = ((v1 << 4) + ((v2 >> 4) & 0xf)) & 0xfff;
      *out16++ = (((v2 & 0xf) << 8 ) + (v3 & 0xff)) & 0xfff;
    }
}

//
// write binary file packing data into 12 bit segments
// this is just an example
// 
void dlardaq::write12(const char *fname, std::vector<adc16_t> &adc)
{
  size_t vecsize = adc.size();
  if(vecsize % 2)
    {
      cout << "The input vector size is not even. Appending zeros" << endl;
      adc.push_back( 0 ); // careful, 0 added returns by reference !!!
      vecsize++;
    }

  // size of the 8 bit buffer
  size_t bufsize = vecsize/2*3;

  // new buffer
  char *buffer = new char[bufsize];

  // fill it
  pack16into12(&adc[0], buffer, bufsize);// in dlardaq.cc
  
  // write it to file
  ofstream fout(fname, ios::out | ios::binary);
  fout.write(buffer, bufsize);
  fout.close();
  
  delete [] buffer;
}


//
// read binary file uppacking 12 bit data into 16 bit segments
//
void dlardaq::read12(const char *fname, std::vector<adc16_t> &adc)
{
  //std::ifstream fin(fname, ios::in | ios::binary | ios::ate);
  std::ifstream fin(fname, ios::in | ios::binary);

  // get number of bytes
  //size_t bufsize = fin.tellg();
  //fin.seekg(0, ios::beg);
  
  // more presise ???
  fin.ignore( std::numeric_limits<std::streamsize>::max() );
  size_t bufsize = fin.gcount();
  fin.clear();  
  fin.seekg( 0, std::ios_base::beg );//c++ func

  size_t npad = 0;
  while( (bufsize+npad) % 3 ) 
    {
      cout << "The input vector size is not divisible by 3. Appending zeros" << endl;
      npad++;
    }

  // allocate memory for file content
  char* buffer = new char[bufsize + npad];
  for(size_t i=bufsize;i<bufsize+npad;i++)
    buffer[i] = '\0';

  // read content of the file
  fin.read(buffer, bufsize);
  fin.close();

  adc.clear();
  adc.resize( (bufsize + npad)/3*2, 0 );
  unpack12into16( buffer, &adc[0], bufsize + npad );// in dlardaq.h

  delete [] buffer;
}


//
// decode run header bytes
ssize_t dlardaq::decode_runhead(const char *buf, dlardaq::runheader_t &rh)
{
  ssize_t rval = RunHeadSz;
  memset(&rh, 0, sizeof(rh));//c++ func
  // in general we will get undefined behaviour if buf < required bytes
  // the try - catch block will not help, but does not hurt either
  try
    {
      // decode run number
      rh.run_num   = (((uint32_t)buf[0] << 24) & 0xFF000000);
      rh.run_num  += (((uint32_t)buf[1] << 16) & 0xFF0000);
      rh.run_num  += (((uint32_t)buf[2] << 8) & 0xFF00);
      rh.run_num  += ( (uint32_t)buf[3] & 0xFF);
      
      rh.run_flags = (uint8_t)buf[4];
    }
  catch(std::exception &e)
    {
      cout << "ERROR: Decode runheader exception " << e.what() << endl;
      rval = -1;
    }
  
  return rval;
}

//
//
// event header bytes
ssize_t dlardaq::decode_evehead(const char* buf, dlardaq::evheader_t &eh)
{
  ssize_t rval = EveHeadSz;
  memset(&eh, 0, sizeof(eh));//c++ func

  try
    {
      // check for delimiting words
      BYTE k0 = buf[0];// in dlardaq.h
      BYTE k1 = buf[1];// in dlardaq.h
      if( !( ((k0 & 0xFF) == EVSKEY) && ((k1 & 0xFF) == EVSKEY) ) )
	{
	  cout << "ERROR: Event delimiting word could not be detected " << endl;
	  throw fex;// in dlardaq.h
	}
      rval = 2;
      // this is actually written in host byte order
      eh.trig_info = ConvertToValue<trigger_t>(buf+rval);
      rval += sizeof(eh.trig_info);
      
      // data quality flags
      eh.dq_flag = (uint8_t)buf[rval++];
      
      // event number 4 bytes
      eh.ev_num   = (((uint32_t)buf[rval++] << 24) & 0xFF000000);
      eh.ev_num  += (((uint32_t)buf[rval++] << 16) & 0xFF0000);
      eh.ev_num  += (((uint32_t)buf[rval++] << 8) & 0xFF00);
      eh.ev_num  += ( (uint32_t)buf[rval++] & 0xFF)-1;


      // event number 4 bytes
      eh.ev_size   = (((uint32_t)buf[rval++] << 24) & 0xFF000000);
      eh.ev_size  += (((uint32_t)buf[rval++] << 16) & 0xFF0000);
      eh.ev_size  += (((uint32_t)buf[rval++] << 8) & 0xFF00);
      eh.ev_size  += ( (uint32_t)buf[rval++] & 0xFF);
    }
  catch(exception &e)
    {
      cout << "ERROR: Decode event header exception "<<e.what() << endl;
      rval = -1;
    }

  return rval;
}

//
// decode run footer to get number of events
ssize_t dlardaq::decode_filefoot(const char *buf, dlardaq::footer_t &rf)
{
  ssize_t rval = FileFootSz;
  memset(&rf, 0, sizeof(rf));//c++ func

  try
    {
      // check for delimiting words
      if( !( ((buf[0]&0xFF) == ENDKEY) && ((buf[1]&0xFF) == ENDKEY) ) )// in dlardaq.h
	{
	  cout << "ERROR: Delimiting word in footer could not be detected " << endl;
	  throw fex;// in dlardaq.h
	}
      rf.num_events = (((uint16_t)buf[2] << 8) & 0xff00) + (buf[3] & 0xff);
    }
  catch(exception &e)//std::exception
    {
      cout << "ERROR: Decode footer info exception "<<e.what() << endl;
      rval = -1;
    }

  return rval;
}



//
// get basic detector information from the decoded run header
ssize_t dlardaq::get_datainfo( const dlardaq::runheader_t &rh, dlardaq::datainfo_t &di )
{
  uint8_t flag = rh.run_flags;
  flag >>= 4; // first 4 bits

  ssize_t rval = -1;

  di.name    = "";
  di.bitadc  = 0;
  di.nch     = 0;
  di.nsample = 0;

  if( flag == 3 ) // 3x1x1 data
    {
      di.bitadc  = dlardaq::BitsADC;
      di.name    = dlardaq::raw311::name;
      di.nch     = dlardaq::raw311::nch;
      di.nsample = dlardaq::raw311::nsa;
      
      rval = 0;
    }
  else if( flag == 6 ) // 6x6x6 data
    {
      di.bitadc  = dlardaq::BitsADC;
      di.name    = dlardaq::raw666::name;
      di.nch     = dlardaq::raw666::nch;
      di.nsample = dlardaq::raw666::nsa;

      rval = 0;
    }
  else
    {
      cout << "ERROR: Uknown detector data type " << (uint16_t)(flag) << endl; 
    }
  
  return rval;
}

