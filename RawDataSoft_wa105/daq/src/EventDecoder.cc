//////////////////////////////////////////////////////////////////////////////////////
//
//  Decoder for the DAQ events for charg readout
//  The general idea is to try to make this object threadsave when possible
//
// 
//////////////////////////////////////////////////////////////////////////////////////

//#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

#include "EventDecoder.h"
#include "LogMsg.h"
#include "Timer.h"

using namespace std;
using namespace dlardaq;

//
// ctor
//
EventDecoder::EventDecoder()
{
  Init();// in EventDecoder
}

//
// ctor -> this method is to be removed, since
// the number of channels and samples should be
// obtained from the run header now
//
EventDecoder::EventDecoder(size_t nch, size_t nsample)
{
  Init();// in EventDecoder
  
  m_di.nch     = nch;// in EventDecoder
  m_di.nsample = nsample;// in EventDecoder
  
  // this is needed to read compressed data
  // m_nadc    = dlardaq::BitsADC;  // bits for ADC resolution
  // m_nch     = nch;                 // total number of channels
  // m_nsample = nsample;             // number of samples in each channel
  
  // m_RunHeadBuf.resize( dlardaq::RunHeadSz );
  // m_EveHeadBuf.resize( dlardaq::EveHeadSz  );
  // m_EndFootBuf.resize( dlardaq::FileFootSz );
  
  // //
  // m_bytes_left = 0;

  // // data mutex initialization
  // pthread_mutex_init(&m_data_mutex, NULL);
}

//
//
//
void EventDecoder::Init()
{
  // this is needed to read compressed data
  //m_nadc    = dlardaq::BitsADC;  // bits for ADC resolution
  //m_nch     = 0;                 // total number of channels
  //m_nsample = 0;                 // number of samples in each channel
  m_di.name    = "";// in EventDecoder
  m_di.bitadc  = dlardaq::BitsADC;  // in EventDecoder. in dlardaq. bits for ADC resolution
  m_di.nch     = 0;// in EventDecoder
  m_di.nsample = 0;// in EventDecoder

  
  m_RunHeadBuf.resize( dlardaq::RunHeadSz );// in EventDecoder. in dlardaq
  m_EveHeadBuf.resize( dlardaq::EveHeadSz  );// in EventDecoder. in dlardaq
  m_EndFootBuf.resize( dlardaq::FileFootSz );// in EventDecoder. in dlardaq
  
  //
  m_bytes_left = 0;// in EventDecoder

  // data mutex initialization
//  pthread_mutex_init(&m_data_mutex, NULL);
}


//
// dtor
//
EventDecoder::~EventDecoder()
{
  Close();// in EventDecoder
//  pthread_mutex_destroy(&m_data_mutex);
}


//
// close input file
//
void EventDecoder::Close()
{
  // lock mutex
//  lock( m_data_mutex );

  if(m_file.is_open())// in EventDecoder
    m_file.close();  // in EventDecoder
  
  m_totev = 0;// in EventDecoder
  
  // unlock mutex
//  unlock( m_data_mutex );
}



//
//
//
ssize_t EventDecoder::Open(std::string finname)
{
  // attempt to close any previously opened files
  Close(); //in EventDecoder
  
  m_EveData.clear(); // in EventDecoder. not used when reading from file

  // clear event bookmarks
  m_events.clear();// in EventDecoder.

  //
  m_file.open(finname.c_str(), ios::in | ios::binary);// in EventDecoder.
  
  if( !m_file.is_open() )// in EventDecoder.
    {
      cout << "Could not open "<<finname<<" for reading" << endl;
//      unlock( m_data_mutex );
      return -1;
    }
  
  // read run header bytes
  ReadBytes( m_RunHeadBuf );// both in EventDecoder.
  dlardaq::decode_runhead(&m_RunHeadBuf[0], m_rnh);// in dlardaq, then both in EventDecoder.
  
  // get data type information
  if( dlardaq::get_datainfo( m_rnh, m_di ) != 0 )// in dlardaq
    {
      cout << "Could not decode detector data type info from "<<finname << endl;
//      unlock( m_data_mutex );
      return -1;
    }

  //
  m_pstart = m_file.tellg();// both in EventDecoder, then c++ func
  
  // read footer
  // fast forward to the end
  m_file.seekg(-dlardaq::FileFootSz, m_file.end );// in EventDecoder, c++ func, in dlardaq, in EventDecoder.
  m_pend = m_file.tellg();// both in EventDecoder., then c++ func
  
  ReadBytes( m_EndFootBuf );//both  in EventDecoder.
  dlardaq::decode_filefoot(&m_EndFootBuf[0], m_flf);// in dlardaq, both in EventDecoder.

  // rewind back to the begining
  m_file.seekg( m_pstart );// in EventDecoder, then c++ func

  
  //
  m_totev = m_flf.num_events;// both in EventDecoder.
  
  if( m_totev == 0 )// in EventDecoder.
    {
      cout << "File "<<finname<<" is empty" << endl;
//      unlock(m_data_mutex);
      Close();// in EventDecoder.
      return 0;
    }

//  unlock( m_data_mutex );
  
  return m_totev;
}

//
// read a given number of bytes from file
//
void EventDecoder::ReadBytes( std::vector<BYTE> &bytes )
{
  m_file.read( &bytes[0], bytes.size() );// in EventDecoder, then c++ func
}


//
//
//
void EventDecoder::ReadEvent( std::vector<adc16_t> &adc, bool headonly )
{
  // first read the header
  ReadBytes( m_EveHeadBuf ); //bot hin EventDecoder
  decode_evehead(&m_EveHeadBuf[0], m_evh);  //in dlardaq, then both in EventDecoder
  //dlardaq::cout << ">>> Event number "<<m_evh.ev_num << endl;
  //dlardaq::cout << ">>> Event size   "<<m_evh.ev_size << endl;
  //dlardaq::cout << ">>> Trig number  "<<m_evh.trig_info.num << endl;
  //dlardaq::cout << ">>> Time stamp   "<<m_evh.trig_info.ts.tv_sec<<" s "
  //<<m_evh.trig_info.ts.tv_nsec<<" ns" << endl;
  //dlardaq::cout << ">>> Flags "<<bitset<8>(m_evh.dq_flag) << endl;

  
  if(headonly) // read only header and skip to next event
    {
      // get event size
      size_t evsz = m_evh.ev_size;// in EventDecoder.

      // move to next event
      m_file.seekg(evsz, std::ios::cur);// in EventDecoder, then c++ func
      return;
    }

  //
  // decode event data
  if(!m_EveDataBuf.empty()) m_EveDataBuf.clear();// in EventDecoder.
  
  m_EveDataBuf.resize(m_evh.ev_size, '\0');//both in EventDecoder.
  
  // read content of the file
  ReadBytes(m_EveDataBuf);// in EventDecoder.
  Decode( &m_EveDataBuf[0], m_EveDataBuf.size(), GETDCFLAG(m_evh.dq_flag), adc );//Decode: in EventDecoder. m_EveDataBuf: in EventDecoder. GETDCFLAG: dlardaq. m_evh: EventDecoder.
}

//
// get a specific event from file
// 
ssize_t EventDecoder::GetEvent(size_t evnum, dlardaq::evheader_t &eh, 
			       std::vector<adc16_t> &adc) 
{
  adc.clear();
  if(!m_file.is_open()) return -1;// in EventDecoder
  if(evnum >= m_totev)  return -1; // in EventDecoder. no-can-do
  
//  lock(m_data_mutex);

  if(evnum < m_events.size())  // in EventDecoder. already know where it is
    {
      m_file.seekg(m_events[evnum]);// in EventDecoder and c++ func
      ReadEvent( adc );// in EventDecoder
      eh = m_evh;
    }
  else
    {
      size_t curev = m_events.size();// in EventDecoder
      while( curev <= evnum )
	{
	  // our event begins at this position
	  streampos pos = m_file.tellg();// streampos: std
	  m_events.push_back( pos );// in EventDecoder
	  
	  // readonly header to get event data size
	  if( curev < evnum ) 
	    ReadEvent(adc, true);// in EventDecoder
	  else // event we want
	    ReadEvent(adc, false);// in EventDecoder
	  
	  curev++;
	}
      
      // last event is the one with want ...
      eh = m_evh;// in EventDecoder
    }

//  unlock(m_data_mutex);

  return evnum;
}

//
//
//
ssize_t EventDecoder::GetEvent( dlardaq::evheader_t &eh,
				std::vector<adc16_t> &adc )
{
  if(m_file.is_open())// in EventDecoder
    {
      // return first event
      return GetEvent( 0, eh, adc);// in EventDecoder
    }

  ssize_t rval = -1;
  //else // we have an event buffer
  if( m_EveData.empty() ) return -1;// in EventDecoder

  // lock mutex
//  lock(m_data_mutex);

  eh  = m_evh;  // in EventDecoder
  adc = m_EveData;// in EventDecoder

  // clear the internal data buffer
  m_EveData.clear();// in EventDecoder
  
  rval = eh.ev_num;
  
  // unlock mutex
//  unlock(m_data_mutex);
  return rval;
}

//
//
//
ssize_t EventDecoder::Decode( const char *buf, size_t nb, bool cflag, 
			      std::vector<adc16_t> &adc )
{
  adc.clear();
  
  //Timer::GetTimer().start();

  if( !cflag ) //not compressed
    {
      adc.resize( (nb*8)/m_di.bitadc, 0 );// in EventDecoder
      unpack12into16( buf, &adc[0], nb );//in dlardaq
    }
  else
    {
      //size_t byteidx = 0;
      //HuffDataCompressor::Instance().DecompressEventData( m_di.bitadc, m_di.nch, 
      //m_di.nsample, buf, nb, byteidx, adc );
      dlardaq::DecompressEventDataFast( buf, nb, m_di.nsample, adc );//in dlardaq, in EventDecoder
      if( adc.size() != m_di.nch * m_di.nsample )// both in EventDecoder
	cout << "The size of the decompressed data does not match the expected" << endl;
    }  

  //Timer::GetTimer().stop();

  return adc.size();
}


//
// we assume that the buffer is 
// run header + event header + event data but it can come in different packets
//
void EventDecoder::ReadBuffer(const char *buf, size_t nb)
{
  if(m_file.is_open())// in EventDecoder
    {
      cout << "Cannot use this function while reading data from a file" << endl;
      return;
    }
  
  // lock mutex
//  lock(m_data_mutex);
  m_totev = 0;  // in EventDecoder
  if(!m_EveData.empty()) m_EveData.clear();// both in EventDecoder
  
  // we have finished reading previous event
  if( m_bytes_left <= 0 || // in EventDecoder
      Timer::GetTimer().splittime(false, false) > TMAXWAIT)// in Timer.h, in EventDecoder
    {
      m_EveDataBuf.clear();// in EventDecoder
      if(IsFirstPacket(buf, nb)) // in EventDecoder. check we have a correct packet
	{
	  // decode header
	  size_t  nb_read = 0;
	  const char *pdata = buf;
	  
	  // decode run header
	  ssize_t ret;
	  ret = dlardaq::decode_runhead(pdata, m_rnh); //in dlardaq
	  if( ret <= 0) 
	    {
	      // unlock mutex
//	      unlock(m_data_mutex);
	      return;
	    }
	  
	  
	  // get detector info
	  ret = dlardaq::get_datainfo( m_rnh, m_di );//in dlardaq
	  if( ret < 0)
	    {
	      // unlock mutex
//	      unlock(m_data_mutex);
	      return;
	    }
	  
	  // move byte inde
	  pdata   += dlardaq::RunHeadSz;//in dlardaq
	  nb_read += dlardaq::RunHeadSz;//in dlardaq
	  //pdata   += ret;
	  //nb_read += ret;
	  

	  // decode event header
	  ret = dlardaq::decode_evehead(pdata, m_evh);//in dlardaq
	  if( ret <= 0) 
	    {
	      // unlock mutex
//	      unlock(m_data_mutex);
	      return;
	    }
	  
	  // move byte index
	  pdata   += dlardaq::EveHeadSz;//in dlardaq
	  nb_read += dlardaq::EveHeadSz;//in dlardaq
	  //pdata   += ret;
	  //nb_read += ret;

	  
	  m_bytes_left = m_evh.ev_size;// both in EventDecoder
	  if(m_bytes_left <= 0) 
	    {
	      // unlock mutex
//	      unlock(m_data_mutex);
	      return;
	    }
	  
	  //cout << "Start recieving event "<<m_evh.ev_num << endl;

	  // copy the data
	  m_EveDataBuf.insert(m_EveDataBuf.end(), pdata, pdata + nb-nb_read);// in EventDecoder
	  m_bytes_left -= (nb - nb_read);// in EventDecoder
	  	  
	  // start timer 
	  Timer::GetTimer().start();// in Timer.h

	  //
//	  unlock(m_data_mutex);
	  return;
	}
      else
	{
//	  unlock(m_data_mutex);
	  return;
	}
    }
  
  // append to buffer
  m_EveDataBuf.insert(m_EveDataBuf.end(), buf, buf+nb);// in EventDecoder
  m_bytes_left -= nb;// in EventDecoder
  if( m_bytes_left < 0 )// in EventDecoder
    {
      cout << "Byte packet mismatch" << endl;
      
      // unlock mutex
//      unlock(m_data_mutex);
      return;
    }

  // we have full event
  if( m_bytes_left == 0 )// in EventDecoder
    {
      // decode data
      bool cflag =  GETDCFLAG(m_evh.dq_flag);//in dlardaq
      Decode( &m_EveDataBuf[0], m_EveDataBuf.size(), cflag, m_EveData );// in EventDecoder
      m_totev = 1;// in EventDecoder
      
      cout << "Recieved: run " << m_rnh.run_num << " event " << m_evh.ev_num << endl;// both in EventDecoder
    }
    
  // unlock mutex
//  unlock(m_data_mutex);
}


//
// detect new event key after run header
//
bool EventDecoder::IsFirstPacket(const char *buf, size_t nb)
{
  BYTE k0 = buf[dlardaq::RunHeadSz];// BYTE: in dlardaq.h.
  BYTE k1 = buf[dlardaq::RunHeadSz+1];
  bool rval = ( ((k0 & 0xFF) == EVSKEY) && ((k1 & 0xFF) == EVSKEY) );// in dlardaq.h.
  
  return rval;
}
