#include "ChMapInterface.h"
#include "LogMsg.h"

using namespace std;
using namespace dlardaq;


//
// ctor
//
ChMapInterface::ChMapInterface()
{
  InitMappings();
}


//
// dtor
//
ChMapInterface::~ChMapInterface()
{
  map<string, dlardaq::ChannelMap*>::iterator it = m_maps.begin();
  for(;it!=m_maps.end();it++)
    delete it->second;
  
  m_maps.clear();
}

//
// add different channel map objects here
//
void ChMapInterface::InitMappings()
{
  // 3x1x1
  m_maps[dlardaq::raw311::name] = new ChannelMap311();
  
  // 6x6x6
  m_maps[dlardaq::raw666::name] = new ChannelMap666();
}


//
//
//
dlardaq::ChannelMap* ChMapInterface::GetMap(std::string key)
{
  map<string, dlardaq::ChannelMap*>::iterator it = m_maps.find(key);
  if( it == m_maps.end() )
    {
      msg_err<<"No channel map for "<<key<<" was found"<<endl;
      return NULL;
    }
  
  return it->second;  
}
