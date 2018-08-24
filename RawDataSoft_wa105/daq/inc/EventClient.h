//////////////////////////////////////////////////////////////////////////////////////////
//
//  Simple event client. It needs to be compiled with
//  CFLAGS -lpthread -levent -levent_pthreads
//  The basic idea is to launch event dispatch in a separate thread so the main thread
//  should be non-blocking. In the case of server disconnect we catch the exception
//  in the signal handler function and try perioducally to reconnect. 
//  To check if the server is alive we periodically send a short message (ping)
//
//  Created: vgalymov@ipnl.inp2p3.fr 
//  Modified: 
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __EVENTSERVER_H__
#define __EVENTSERVER_H__

// only included with WITHLIBEVENT2 preprocessor is set
#ifdef WITHLIBEVENT2

#include <map>
#include <stdint.h>
#include <pthread.h>
#include <arpa/inet.h>
#include <event2/event.h>
#include <event2/thread.h>
#include <event2/buffer.h>
#include <event2/bufferevent.h>

#include "EventDecoder.h"

#define RWTIMEOUT_SEC 10   // r/w timeout in seconds
#define RCTIMEOUT_SEC 15   // ping / reconnect timeout
//
namespace dlardaq
{
  class EventClient
  {
  public:
    EventClient();
    ~EventClient();
    
    // connect to data server
    void Connect(uint16_t port, const char *serv_ip);

    // check if we are up
    bool Connected() const { return m_connected; }

    //
    void SetDecoder(dlardaq::EventDecoder *ptr){ m_decoder = ptr; }

    // set a pointer to the event buffer
    //void SetEventBuffer();

    // stop connection thread
    void Stop();
    
  private:
    // thread function
    static void* threadProcess(void *obj)
    {
      ((EventClient *)obj)->runThread();
      return NULL;
    }

    void setConnected(bool val){ m_connected = val; }
    
    // !!! DO NOT call this function in callbacks as it !!!
    // !!! will delete event_base and cause segfault    !!! 
    void Cleanup();
    
    // try to send dummpy message to server
    // to check if it is still alive  connected
    void ping();

    //
    void runThread();
    void stopThread();

    // timer add / remove
    //void addTimer(){ event_add(ev_timer, &rc_timer); }
    //void removeTimer(){ event_del(ev_timer); }

    // add remove buffer event
    void addBufEvent();
    void removeBufEvent();

    // read event data
    void readEvent();

    // check connection status periodically
    static void connect_cb(evutil_socket_t fd, short what, void *arg);
    
    // callback for signals (sigpipe)
    static void signal_cb(evutil_socket_t sig, short what, void *arg); 
    
    // buffer callback functions
    static void buffer_readcb(struct bufferevent *bev, void *arg);
    static void buffer_writecb(struct bufferevent *bev, void *arg);
    static void buffer_eventcb(struct bufferevent *bev, short what, void *arg);

    // decoder for event data
    EventDecoder* m_decoder;

    // flags
    bool m_running;
    bool m_connected;
    
    // server net stuff
    struct sockaddr_in m_sin;
    
    // run thread
    pthread_t m_run_thread;
        
    // lib event
    struct timeval      rc_timer;
    struct event*       ev_timer;
    struct event_base*  ev_base;
    struct bufferevent* bev;
  };  
};

#endif //WITHLIBEVENT2

#endif //ifdef



