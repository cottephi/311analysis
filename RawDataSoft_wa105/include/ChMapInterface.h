////////////////////////////////////////////////////////
//
//  this is singleton class
//  Add different ch mappings to InitMappings() function
//
//  
//
////////////////////////////////////////////////////////



#ifndef __CHMAPINTEFACE_H__
#define __CHMAPINTEFACE_H__

#include <map>
#include <string>
#include "ChannelMap.h"

namespace dlardaq
{
  class ChMapInterface
  {
  public:
    static ChMapInterface& Instance()
    {
      static ChMapInterface inst;
      return inst;
    }
    
    ChannelMap* GetMap(std::string key);
    
  private:
    // ctor
    ChMapInterface();
    ChMapInterface(const ChMapInterface&);
    ChMapInterface& operator=(const ChMapInterface&);
    // Prevent unwanted destruction
    ~ChMapInterface();

    void InitMappings();
    
    // channel maps for different geometries
    std::map<std::string, ChannelMap*> m_maps;
    
  };  
};


#endif
