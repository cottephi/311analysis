/////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//  Mapping of the DAQ channels to CRP views
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __CHANNELMAP_H__
#define __CHANNELMAP_H__

#include <string>
#include <vector>

#include "dlardaq.h"

namespace dlardaq
{

  //
  //
  //
  typedef struct daqchan_t
  {
    int crate;
    int card;
    int ch;
  } daqchan_t;

  typedef std::vector<daqchan_t> daqchans;

  //
  // channels in a given CRP module
  //
  class ViewChansCRM
  {
  public:
    ViewChansCRM(){;}
    ViewChansCRM(std::vector<size_t> &nch_view);
    ~ViewChansCRM(){;}
  
    int Set(int iview, int chv, daqchan_t &daqch);
    int Get(int iview, int chv, daqchan_t &daqch);
    int IsSet( int iview, int chv )
    {
      if(!CheckSize(iview, chv)) return -1;
      if(m_views[iview][chv].crate<0) return 0;
      return 1;   
    }


  private:
    bool CheckSize( int iview, int chv )
    {
      if( iview < 0 || iview >= (int)m_views.size() ) return false;
      if( chv < 0 || chv >= (int)m_views[iview].size() ) return false;
      return true;
    }

    std::vector<daqchans> m_views;
  };

  //
  //
  //
  class ChannelMap
  {
  public:
  ChannelMap(): 
    m_MaxCRMs(0), m_MaxViews(0), m_MaxChPerView(0),
      m_Crates(0), m_CardPerCrate(0), m_ChPerCard(0){;}
    virtual ~ChannelMap(){;}
  
    virtual int MapToCRP(int crate, int card, int ch, int &crm, int &view, int &chv) = 0;
    virtual int MapToCRP(int seqch, int &crm, int &view, int &chv) = 0;
    virtual int MapToDAQ(int crm, int view, int chv, int &crate, int &card, int &ch) = 0;
    virtual int MapToDAQFile(int crm, int view, int chv) = 0;
    //virtual int WriteAscii(const char *fname) = 0;
  
  protected:
    virtual bool CheckDaqInputs(int crate, int card, int ch);
    virtual bool CheckCRPInputs(int crm, int view, int ch);
  
    // CRP parameters
    size_t m_MaxCRMs;
    size_t m_MaxViews;
    size_t m_MaxChPerView;

    // DAQ parameters
    size_t m_Crates;
    size_t m_CardPerCrate;
    size_t m_ChPerCard;
    size_t m_ChPerCrate; 
  };

  //
  // for 3x1x1
  //
  class ChannelMap311 : public ChannelMap
  { 
  public:
    ChannelMap311();
    ~ChannelMap311(){;}
  
    int MapToCRP(int crate, int card, int ch, int &crm, int &view, int &chv);
    int MapToCRP(int seqch, int &crm, int &view, int &chv);
    int MapToDAQ(int crm, int view, int chv, int &crate, int &card, int &ch);
    // this function should return a ch number in file sequence
    int MapToDAQFile(int crm, int view, int ch); 
    //int WriteAscii(const char *fname);

  private:
    void InitMap();
  
    size_t m_Half;
    size_t m_Offset; 
    size_t m_ChPerView0;
    size_t m_ChPerView1;

    ViewChansCRM m_VwChans;
  };


  //
  // for 6x6x6
  //
  class ChannelMap666 : public ChannelMap
  { 
  public:
    ChannelMap666();
    ~ChannelMap666(){;}
  
    int MapToCRP(int crate, int card, int ch, int &crm, int &view, int &chv);
    int MapToCRP(int seqch, int &crm, int &view, int &chv);
    int MapToDAQ(int crm, int view, int chv, int &crate, int &card, int &ch);
    // this function should return a ch number in file sequence
    int MapToDAQFile(int crm, int view, int ch); 
    //int WriteAscii(const char *fname);

  private:
    void InitMap();
  
    size_t m_ChPerView0;
    size_t m_ChPerView1;

    //
    size_t m_HalfCrate;
    size_t m_QuartCrate;
    size_t m_HalfCard;
    
    std::vector< ViewChansCRM > m_VwChans;
  };
};

#endif
