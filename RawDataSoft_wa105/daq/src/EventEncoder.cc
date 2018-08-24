////////////////////////////////////////////////////////////////////////
//
//
//  Class to produce raw data files / events  for testing purposes
//
//
////////////////////////////////////////////////////////////////////////

//#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string.h>

#include "EventEncoder.h"
#include "LogMsg.h"

using namespace std;
using namespace dlardaq;


// 
// ctor
//
EventEncoder::EventEncoder(size_t nch, size_t nsample)
{
  m_nadc    = dlardaq::BitsADC; // bits for ADC
  m_nch     = nch;              // number of channels
  m_nsample = nsample;          // number of samples in each channel
  
  // expected data size in bytes
  m_datasz  = m_nadc * m_nch * m_nsample; //data size in bits
  if(m_datasz % 2) m_datasz += 8;  // TODO: warning
  m_datasz /= 8; //in bytes

  m_RunHeadBytes = dlardaq::RunHeadSz;
  m_EveHeadBytes = dlardaq::EveHeadSz;
  m_EndFootBytes = dlardaq::FileFootSz;

  m_RunHeadBuf = new BYTE[m_RunHeadBytes];
  m_EveHeadBuf = new BYTE[m_EveHeadBytes];
  m_EndFootBuf = new BYTE[m_EndFootBytes];
  
  m_totev     = 0;
  m_totbytes  = 0;
  m_runnb     = 6102016;
  m_dcflag    = false;
}


//
// dtor
//
EventEncoder::~EventEncoder()
{
  Close();

  delete [] m_RunHeadBuf;
  delete [] m_EveHeadBuf;
  delete [] m_EndFootBuf;
}


//
// close input file
//
void EventEncoder::Close()
{
  if(fOut.is_open())
    {
      // write out total number of events
      EncodeRunFoot();
      WriteBytes( m_EndFootBuf, m_EndFootBytes );
      fOut.close();  
    }
  
  m_totev = 0;
  m_totbytes  = 0;
}


//
// open file for writing
//
int EventEncoder::Open(std::string finname, uint16_t flag)
{
  // attempt to close any previously opened files
  Close();
  
  m_totev = 0;
  m_totbytes  = 0;

  //
  fOut.open(finname.c_str(), ios::out | ios::binary);
  
  if( !fOut.is_open() )
    {
      msg_err<<"Could not open "<<finname<<" for writing"<<endl;
      return -1;
    }

  // write run header bytes
  EncodeRunHead( flag );
  WriteBytes( m_RunHeadBuf, m_RunHeadBytes );
  
  return 1;
}

//
//
//
void EventEncoder::WriteBytes(const BYTE* buf, size_t nb)
{
  m_totbytes += nb;

  fOut.write(buf, nb); 
}

//
//
//
void EventEncoder::WriteEvent(std::vector<adc16_t> &adc)
{
  // make event header
  vector<BYTE> buf;  
  EncodeEvent(adc, buf);
  WriteBytes( &buf[0], buf.size() );
  buf.clear();
}

//
//
//
void EventEncoder::WriteEvent(float ped, float rms)
{
  std::vector<adc16_t> adc;
  GenerateRandom( ped, rms, adc );
  WriteEvent( adc );
  adc.clear();
}

//
//
//
void EventEncoder::WriteBufferEvent(float ped, float rms, std::vector<BYTE> &bytes)
{
  std::vector<adc16_t> adc;
  GenerateRandom( ped, rms, adc );
  WriteBufferEvent( adc, bytes );
  adc.clear();
}

//
//
//
void EventEncoder::WriteBufferEvent(std::vector<adc16_t> &adc, std::vector<BYTE> &bytes)
{
  bytes.clear();
  EncodeRunHead( 0x30 ); //
  bytes.insert(bytes.begin(), m_RunHeadBuf, m_RunHeadBuf+m_RunHeadBytes);
  EncodeEvent( adc, bytes );
}

//
//
//
void EventEncoder::EncodeEvent( std::vector<adc16_t> &adc, 
				std::vector<BYTE> &buf )
{
  uint32_t evsize = 0;

  size_t start = buf.size();
  buf.resize( start + m_EveHeadBytes );

  if(!m_dcflag)
    {
      size_t sz = buf.size() + m_datasz; 
      buf.resize( sz, '\0' );
      
      // fill it
      pack16into12(&adc[0], &buf[start+m_EveHeadBytes], m_datasz);
      evsize = m_datasz;
      m_totbytes += m_datasz;
    }
  else // compress the event data
    {
      vector<BYTE> tmp;
      HuffDataCompressor::Instance().CompressEventData(m_nadc, m_nch, m_nsample, adc, tmp);
      evsize = tmp.size();
      buf.insert(buf.end(), tmp.begin(), tmp.end());
    }

  EncodeEventHead(evsize);
  for(size_t i=0;i<m_EveHeadBytes;i++)
    buf[start+i] = m_EveHeadBuf[i];
  
  m_totev++;
}

//
// dummy run header
//
void EventEncoder::EncodeRunHead( uint16_t rflag )
{

  size_t pos = 0;
  m_RunHeadBuf[pos++] = (m_runnb >> 24) & 0xff;
  m_RunHeadBuf[pos++] = (m_runnb >> 16) & 0xff;
  m_RunHeadBuf[pos++] = (m_runnb >> 8) & 0xff;
  m_RunHeadBuf[pos++] = (m_runnb) & 0xff;
  m_RunHeadBuf[pos++] = rflag & 0xff;
  //m_RunHeadBuf[pos++] = 0x11;
}


//
// dummy event header
//
void EventEncoder::EncodeEventHead(uint32_t evsize)
{
  size_t pos = 0;
  m_EveHeadBuf[pos++] = EVSKEY;
  m_EveHeadBuf[pos++] = EVSKEY;


  trigger_t trig;
  memset(&trig, 0, sizeof(trig));
  
  trig.ts.tv_sec  = time(NULL);
  trig.ts.tv_nsec = rand() % 125000000;
  trig.type       = rand() % 3 + 1;
  trig.num        = m_totev;
  
  memcpy(m_EveHeadBuf+pos,(void*)&trig, sizeof(trig));
  pos += sizeof(trig);
  uint8_t dqflag = 0x0;
  if(m_dcflag) SETDCFLAG(dqflag); // set compression flag
  m_EveHeadBuf[pos++] = dqflag;
  
  // last 4 bytes for event number
  m_EveHeadBuf[pos++] = (m_totev >> 24) & 0xff;
  m_EveHeadBuf[pos++] = (m_totev >> 16) & 0xff;
  m_EveHeadBuf[pos++] = (m_totev >> 8) & 0xff;
  m_EveHeadBuf[pos++] = (m_totev) & 0xff;

  m_EveHeadBuf[pos++] = (evsize >> 24) & 0xff;
  m_EveHeadBuf[pos++] = (evsize >> 16) & 0xff;
  m_EveHeadBuf[pos++] = (evsize >> 8) & 0xff;
  m_EveHeadBuf[pos++] = (evsize) & 0xff;
}


//
// total number of events for footer
//
void EventEncoder::EncodeRunFoot()
{
  size_t pos = 0;
  m_EndFootBuf[pos++] = ENDKEY;
  m_EndFootBuf[pos++] = ENDKEY;
  m_EndFootBuf[pos++] = (m_totev >> 8) & 0xff;
  m_EndFootBuf[pos++] = (m_totev) & 0xff;
}


//
//
//
void EventEncoder::GenerateRandom(float ped, float rms, 
				  std::vector<adc16_t> &adc)
{
  adc.clear();
  adc.resize( m_nch * m_nsample );

  for(size_t i=0;i<m_nch;i++)
    {
      for(size_t j=0;j<m_nsample;j++)
	{
	  if( i==0 )
	    {
	      double val = genNoise( ped, rms );
	      val = floor(val + 0.5);
	      if(val < 0 ) val = 0;
	      adc[j] = (adc16_t)(val);
	    }
	  else
	    {
	      adc[j + i*m_nsample] = adc[j];
	    }
	}
    }
}

//
// Box-Muller transform 
// https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
//
double EventEncoder::genNoise(double mu, double sigma)
{
  const double epsilon = std::numeric_limits<double>::min();
  const double two_pi = 2.0*3.14159265358979323846;
  
  static double z0, z1;
  static bool generate;
  generate = !generate;
  
  if (!generate)
    return z1 * sigma + mu;
  
  double u1, u2;
  do
    {
      u1 = rand() * (1.0 / RAND_MAX);
      u2 = rand() * (1.0 / RAND_MAX);
    }
  while ( u1 <= epsilon );
  
  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}
