/////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//  Mapping of the DAQ channels to CRP views
//  
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <cmath>

#include "ChannelMap.h"
#include "LogMsg.h"

using namespace std;
using namespace dlardaq;
//
//
//
ViewChansCRM::ViewChansCRM(std::vector<size_t> &nch_view)
{
  m_views.resize( nch_view.size() );
  daqchan_t tmp;
  tmp.crate = -1;
  tmp.card  = -1;
  tmp.ch    = -1;
  
  for(size_t i=0;i<m_views.size();i++)
    {
      m_views[i].resize( nch_view[i], tmp );
    }
}


//
//
//
int ViewChansCRM::Set( int iview, int chv, daqchan_t &daqch )
{
  if(!CheckSize( iview, chv ) ) return -1;
  m_views[iview][chv] = daqch;
  return 1;
}

//
//
//
int ViewChansCRM::Get(int iview, int chv,  daqchan_t &daqch)
{
  if(IsSet(iview, chv) <= 0) return -1;
  daqch = m_views[iview][chv];
  return 1;
}


//
//
//
//
bool ChannelMap::CheckDaqInputs(int crate, int card, int ch)
{
  if( crate >= (int)m_Crates || crate < 0)
    {
      dlardaq::msg_err<<"Crate id "<<crate<<" is not correct. Max value "<<m_Crates-1<<endl;
      return false;
    }

  if( card >= (int)m_CardPerCrate || card < 0)
    {
      dlardaq::msg_err<<"Card id "<<card<<" is not correct. Max value "<<m_CardPerCrate-1<<endl;
      return false;
    }

  if( ch >= (int)m_ChPerCard || ch<0 )
    {
      dlardaq::msg_err<<"Ch id "<<ch<<" is not correct. Max value "<<m_ChPerCard-1<<endl;
      return false;
    }

  return true;
}

//
//
//
bool ChannelMap::CheckCRPInputs(int crm, int view, int ch)
{
  if( crm < 0 || crm >= (int)m_MaxCRMs )
    {
      dlardaq::msg_err<<"CRM id "<<crm<<" is not valid"<<endl;
      return false;      
    }
  
  if( view < 0 || view >= (int)m_MaxViews ) 
    {
      dlardaq::msg_err<<"View number "<<view<<" is not valid"<<endl;
      return false;
    }
  
  if( view < 0 || view >= (int)m_MaxChPerView )
    {
      dlardaq::msg_err<<"View channel number "<<ch<<" is not valid"<<endl;
      return false;
    }

  return true;
}

//
//
//
ChannelMap311::ChannelMap311()
{
  //m_Crates       = 4;
  //m_CardPerCrate = 5; 
  //m_ChPerCard    = 64; 

  m_Crates       = dlardaq::raw311::ncrate;
  m_CardPerCrate = dlardaq::raw311::nccrate; 
  m_ChPerCard    = dlardaq::raw311::nchcard; 

  m_ChPerCrate = m_ChPerCard*m_CardPerCrate;
  m_Half       = m_ChPerCard/2;
  m_Offset     = m_CardPerCrate*m_Half;

  m_MaxCRMs  = 1;
  m_MaxViews = 2;
  m_MaxChPerView = 960;
  
  //
  m_ChPerView0 = 320;
  m_ChPerView1 = 960;

  std::vector<size_t> tmp(2);
  tmp[0] = m_ChPerView0;
  tmp[1] = m_ChPerView1;
  m_VwChans = ViewChansCRM(tmp);
  
  InitMap();
}


//
//
//
void ChannelMap311::InitMap()
{
  // view 0;
  int crm, iview, chv;
  size_t count0 = 0;
  daqchan_t dqch;
  dqch.crate = 0;
  for(size_t card=0;card<m_CardPerCrate;card++)
    {
      dqch.card = card;
      for(size_t c=0;c<m_ChPerCard;c++)
	{
	  dqch.ch = c;
	  if(MapToCRP(0, card, c, crm, iview, chv)>0)
	    {
	      //dlardaq::msg_info<<iview<<" "<<chv<<endl;
	      if(m_VwChans.Set(iview, chv, dqch)>0)
		count0++;
	      //else
	      //dlardaq::msg_warn<<iview<<" "<<chv<<endl;
	    }
	}
    }
  
  if(count0!=m_ChPerView0)
    {
      dlardaq::msg_err<<"Incorrect map strucutre in view 0"<<endl;
      exit(-1);
    }


  // view 1
  size_t count1 = 0;
  for(size_t crate=1;crate<m_Crates;crate++)
    {
      dqch.crate = crate;
      for(size_t card=0;card<m_CardPerCrate;card++)
	{
	  dqch.card = card;
	  for(size_t c=0;c<m_ChPerCard;c++)
	    {
	      dqch.ch = c;
	      if(MapToCRP(crate, card, c, crm, iview, chv)>0)
		{
		  if(m_VwChans.Set(iview, chv, dqch)>0)
		    count1++;
		}
	    }
	}
    }

  if(count1!=m_ChPerView1)
    {
      dlardaq::msg_err<<"Incorrect map strucutre in view 1"<<endl;
      exit(-1);
    }
}
  

//
// 3x1x1 channel mappings
//
int ChannelMap311::MapToCRP(int crate, int card, int ch, int &crm, int &view, int &chv)
{
  if( !CheckDaqInputs(crate, card, ch) )
    {
      dlardaq::msg_info<<crate<<", "<<card<<", "<<ch<<std::endl;
      return -1;
    }

  // fix for 8 ch inversion
  int tmp0 = (int)(ch/8);
  ch = (-(ch % 8) + 7) + tmp0 * 8;

  // only 1 crm
  crm = 0;
  
  // views
  view = (crate==0)?0:1;
  
  //
  int tmp1, tmp2, tmp3;
  if( crate == 0) 
    {
      // simple mapping 
      //tmp1 = 0;
      //tmp2 = card * m_ChPerCard;
      //tmp3 = ch;
      
      ch = (-(ch%32) + 31) + (int)(ch/32)*32;
      
      tmp1 = 0;
      tmp2 = ch<(int)m_Half ? 0:m_Offset;
      tmp3 = ch<(int)m_Half ? ch:(ch-m_Half);
    }
  else 
    {
      // simple mapping 
      //tmp1 = (crate-1)*m_ChPerCrate;
      //tmp2 = card * m_ChPerCard;
      //tmp3 = ch;


      tmp1 = (crate - 1)*m_ChPerCrate;
      tmp2 = ch<(int)m_Half ? m_Offset:0;
      tmp3 = ch<(int)m_Half ? ch:(ch-m_Half);
    }
  
  chv = card*m_Half + tmp1 + tmp2 + tmp3; 
  
  // int tmp4 = -(chv % 8) + 7;
  // if( tmp4 < 0 ) 
  //   {
  //     dlardaq::msg_err<<"I screwed up"<<endl;
  //     exit(-1);
  //   }
  // int tmp5 = (int)((float)(chv)/8.0);
  // chv = tmp5 + tmp4;
  
  // !!! make right handed coordinate system !!!
  if( crate == 0 ) chv = -chv + (m_ChPerCrate - 1);
  
  // simple mapping 
  //chv = tmp1 + tmp2 + tmp3; 
  
  return 1;
}

//
// map continuous sequence of daq channels to CRP view channels
//
int ChannelMap311::MapToCRP(int seqch, int &crm, int &view, int &chv)
{
  if( seqch < 0 || seqch >= (int)(m_Crates * m_ChPerCrate) ) return -1;

  int crate = seqch / m_ChPerCrate;
  int card  = (seqch - crate * m_ChPerCrate) / m_ChPerCard;
  int ch    = (seqch - crate * m_ChPerCrate - card * m_ChPerCard);
  
  return MapToCRP( crate, card, ch, crm, view, chv );
}


//
// 
//
int ChannelMap311::MapToDAQ(int crm, int view, int chv, int &crate, int &card, int &ch)
{
  crate = -1;
  card  = -1;
  ch    = -1;
  
  if(!CheckCRPInputs(crm, view, chv)) return -1;

  daqchan_t dqch;
  if(m_VwChans.Get(view, chv, dqch)<0) return -1;
  
  crate = dqch.crate;
  card  = dqch.card;
  ch    = dqch.ch;  

  return 1;
}


//
// return 1d map
//
int ChannelMap311::MapToDAQFile(int crm, int view, int chv)
{
  int crate, card, ch;
  if(MapToDAQ(crm, view, chv, crate, card, ch)<0) return -1;
  
  //
  int seqch = crate * m_ChPerCrate + card*m_ChPerCard + ch;
  return seqch;
}




//
//
//
ChannelMap666::ChannelMap666()
{
  m_Crates       = dlardaq::raw666::ncrate;
  m_CardPerCrate = dlardaq::raw666::nccrate; 
  m_ChPerCard    = dlardaq::raw666::nchcard; 
  m_ChPerCrate   = m_CardPerCrate * m_ChPerCard;

  m_MaxCRMs      = 4;   // number of 3x3 m2 CRMs
  m_MaxViews     = 2;   // number of readout views
  m_MaxChPerView = 960; // number of channels in view


  // channels per view in each 3x3 m2 CRM
  m_ChPerView0 = 960;
  m_ChPerView1 = 960;

  //
  m_HalfCrate  = m_Crates / 2;
  m_QuartCrate = m_HalfCrate / 2;

  //
  m_HalfCard   = m_CardPerCrate / 2;
  
  //
  std::vector<size_t> tmp(2);
  tmp[0] = m_ChPerView0;
  tmp[1] = m_ChPerView1;
  for( size_t i=0;i<m_MaxCRMs;i++)
    m_VwChans.push_back( ViewChansCRM(tmp) );
  
  InitMap();
}



//
//
//
void ChannelMap666::InitMap()
{
  int crm, iview, chv;
  daqchan_t dqch;
  
  size_t count = 0;
  for(size_t crate = 0; crate < m_Crates; ++crate)
    {
      dqch.crate = crate;
      for(size_t card = 0; card < m_CardPerCrate; ++card)
	{
	  dqch.card = card;
	  for(size_t c = 0; c < m_ChPerCard; ++c)
	    {
	      dqch.ch = c;
	      if(MapToCRP(crate, card, c, crm, iview, chv) > 0 )
		{
		  if(crm < (int)m_MaxCRMs)
		    {
		      if( m_VwChans[crm].Set(iview, chv, dqch ) > 0)
			count++;
		    }

		}
	    } // loop ch in card
	} // loop cards in crate
    } // loop crates
  
  if(count != dlardaq::raw666::nch)
    {
      dlardaq::msg_err<<"Incorrect map strucutre for channel mapping of 666 detected"<<endl;
      exit(-1);
    }
}
  

//
// 6x6x6 channel mappings
//
int ChannelMap666::MapToCRP(int crate, int card, int ch, int &crm, int &view, int &chv)
{
  if( !CheckDaqInputs(crate, card, ch) )
    {
      dlardaq::msg_info<<crate<<", "<<card<<", "<<ch<<std::endl;
      return -1;
    }

  // Find the view
  // crate 0 - 5 are view 0 (or X) and 6 - 11 are view 1 (or Y)
  view = (crate<(int)m_HalfCrate)?0:1;

  // Find the CRM number
  // cards 0 - 4 read 3x1 of one CRM and 5 - 9 read the 3x1 of the other
  bool tmp1 = (card  >= (int)m_HalfCard);
  bool tmp2 = ((crate - view * m_HalfCrate) >= m_QuartCrate);
  crm = 0;
  if( view == 0 )
    {
      if( tmp1 ) crm |= 1 << 1;
      if( tmp2 ) crm |= 1;
    }
  else
    {
      if( tmp2 ) crm |= 1 << 1;
      if( tmp1 ) crm |= 1;
    }
  
  
  // Find the view channel
  chv = ch + ( (crate % m_QuartCrate) * m_HalfCard + ( card % m_HalfCard ) ) * m_ChPerCard;
  
  return 1;
}

//
// map continuous sequence of daq channels to CRP view channels
//
int ChannelMap666::MapToCRP(int seqch, int &crm, int &view, int &chv)
{
  if( seqch < 0 || seqch >= (int)(dlardaq::raw666::nch) ) return -1;

  int crate = seqch / m_ChPerCrate;
  int card  = (seqch - crate * m_ChPerCrate) / m_ChPerCard;
  int ch    = (seqch - crate * m_ChPerCrate - card * m_ChPerCard);
  
  return MapToCRP( crate, card, ch, crm, view, chv );
}


//
// 
//
int ChannelMap666::MapToDAQ(int crm, int view, int chv, int &crate, int &card, int &ch)
{
  crate = -1;
  card  = -1;
  ch    = -1;
  
  if(!CheckCRPInputs(crm, view, chv)) return -1;

  daqchan_t dqch;
  if(m_VwChans[crm].Get(view, chv, dqch)<0) return -1;
  
  crate = dqch.crate;
  card  = dqch.card;
  ch    = dqch.ch;  

  return 1;
}


//
// return 1d map
//
int ChannelMap666::MapToDAQFile(int crm, int view, int chv)
{
  int crate, card, ch;
  if(MapToDAQ(crm, view, chv, crate, card, ch)<0) return -1;
  
  //
  int seqch = crate * m_ChPerCrate + card*m_ChPerCard + ch;
  return seqch;
}
