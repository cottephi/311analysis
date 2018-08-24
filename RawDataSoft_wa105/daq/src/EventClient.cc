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

#ifdef WITHLIBEVENT2

//
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/socket.h>

#include "EventClient.h"
#include "LogMsg.h"

using namespace std;
using namespace dlardaq;

//
// ctor
//
EventClient::EventClient()
{
  // enable treading support
  evthread_use_pthreads();
  
  memset(&m_sin, 0, sizeof(m_sin));

  m_running   = false;
  m_connected = false;

  //
  ev_timer = NULL;
  ev_base  = NULL;
  bev      = NULL;

  // pointer to event decoder (if any)
  m_decoder = NULL;
}

//
// dtor
//
EventClient::~EventClient()
{
  Stop();
  Cleanup();
}

//
// connect
//
void EventClient::Connect(uint16_t port, const char *serv_ip)
{
  msg_info<<"Connecting to "<<serv_ip<<":"<<port<<endl;
  if( m_connected && ntohs(m_sin.sin_port) == port )
    {
      msg_warn<<"already connected on port "<<port<<endl;
      return;
    }
  else if(m_running) Stop();
  Cleanup();
  
  m_running   = false;
  m_connected = false;  
  
  m_sin.sin_family = AF_INET;
  m_sin.sin_port   = htons(port);
  if(serv_ip) inet_aton(serv_ip, &m_sin.sin_addr);
  else m_sin.sin_addr.s_addr = htonl(0x7f000001); // localhost
  
  // create event base
  if( (ev_base = event_base_new()) == NULL )
    {
      //err(1,"failed to create event base");
      msg_err<<"failed to create event base"<<endl;
      exit(1);
    }
  

  // timer event
  ev_timer = event_new( ev_base, -1, EV_PERSIST | EV_TIMEOUT, EventClient::connect_cb, this );
  rc_timer.tv_sec = RCTIMEOUT_SEC; rc_timer.tv_usec = 0;
  event_add(ev_timer, &rc_timer);
  
  // buffer event
  m_connected = true;
  bev = bufferevent_socket_new( ev_base, -1, BEV_OPT_CLOSE_ON_FREE | BEV_OPT_THREADSAFE );
  bufferevent_setcb(bev, buffer_readcb, buffer_writecb, buffer_eventcb, this);
  bufferevent_enable(bev, EV_READ | EV_WRITE); 
  if(bufferevent_socket_connect(bev,(struct sockaddr *)&m_sin, sizeof(m_sin)) != 0)
    {
      msg_err<<"failed buffer socket connect"<<endl;
      Cleanup();
      exit(1);
    }
  
  // handle interrupt signal events
  event_add(evsignal_new(ev_base, SIGPIPE, EventClient::signal_cb, this), NULL);
  event_add(evsignal_new(ev_base, SIGINT, EventClient::signal_cb, this), NULL);
  //event_add(evsignal_new(ev_base, SIGHUP, EventClient::signal_cb, this), NULL);

  // start listener thread
  if( pthread_create(&m_run_thread, NULL, threadProcess, this) != 0 )
    {
      msg_err<<"failed to create thread connect"<<endl;
      Cleanup();
      exit(1);
    }

  //
  m_running = true;
}


//
// clean up
//
void EventClient::Cleanup()
{
  memset(&m_sin, 0, sizeof(m_sin));

  if(bev)
    {
      bufferevent_free( bev );
      bev = NULL;
    }
  if(ev_timer)
    {
      event_free( ev_timer );
      ev_timer = NULL;      
    }
  if(ev_base)
    {
      event_base_free( ev_base );
      ev_base = NULL;
    }  
}

//
// this function is non-blocking
//
void EventClient::Stop()
{
  if(!m_running) return;
  
  m_running   = false;
  m_connected = false;
  
  //
  stopThread();
  
  // join the dispatch thread
  pthread_join( m_run_thread, NULL);
}

//
// dispatch event base loop 
//
void EventClient::runThread()
{
  //
  event_base_dispatch(ev_base);
}

//
// send exit to dispatch loop
//
void EventClient::stopThread()
{
  //
  event_base_loopexit(ev_base, NULL);  
}

//
// add buffer event (connect to server)
// TODO: the error handling could be improved here
//
void EventClient::addBufEvent()
{
  bev = bufferevent_socket_new( ev_base, -1, BEV_OPT_CLOSE_ON_FREE | BEV_OPT_THREADSAFE );
  bufferevent_setcb(bev, buffer_readcb, buffer_writecb, buffer_eventcb, this);
  bufferevent_enable(bev, EV_READ | EV_WRITE); 
  bufferevent_socket_connect(bev,(struct sockaddr *)&m_sin, sizeof(m_sin));
}

//
// remove buffer event (e.g., server has diconnected)
//
void EventClient::removeBufEvent()
{
  bufferevent_free( bev );
  bev = NULL;
}

//
// Send a short message to the server periodically
// if it stopped we will get an exception/error 
// which we should catch, close the socket, 
// and try periodically to reconnect
//
void EventClient::ping()
{
  char buf[8] = "hello?";
  bufferevent_write(bev, buf, sizeof(buf));
}

//
// read event data
//
void EventClient::readEvent()
{
  struct evbuffer *input = bufferevent_get_input(bev);
  //size_t bytes_toread   = evbuffer_get_length(input);
  //cout<<"bytes to read "<<bytes_toread<<endl;
  //evbuffer_remove(input, m_buffer, bytes_toread);
  
  ssize_t nb;
  char buf[4096];
  while( (nb = evbuffer_remove(input, buf, sizeof(buf))) > 0 )
    {
      if( m_decoder ) 
	m_decoder->ReadBuffer(buf, nb);
      
    }
}
 


//
// accept_conn_cb
//
void EventClient::connect_cb(evutil_socket_t fd, short what, void *arg)
{
  EventClient* obj = static_cast<EventClient*>(arg);
  if(obj->Connected())  // ping the connection
    {
      //cout<<"Send ping"<<endl;
      obj->ping();
    }
  else
    {
      //cout<<"Remake new"<<endl;
      obj->addBufEvent();
    }
}


//
// signal handling callback
//
void EventClient::signal_cb(evutil_socket_t sig, short what, void *arg)
{
  
  // signal handling of different interrupts
  
  switch( sig )
    {
    case SIGPIPE: 
      //warn("sigpipe");
      // ignore broken pipe and and let buffer_eventcb
      // take care of closing disconnected server
      break;
      
    case SIGINT:
      (static_cast<EventClient*>(arg))->Stop();
      break;
      
    default:
      //warn("No default action on interrupt signal %d", (int)sig);
      msg_warn<<"No default action on interrupt signal "<<sig<<endl;
      break;      
    }
}


//
// drain the input buffer
//
void EventClient::buffer_readcb(struct bufferevent *bev, void *arg)
{
  (static_cast<EventClient*>(arg))->readEvent();
}


//
// evbuffer_get_length(bufferevent_get_output(bev)) should be 0
// 
void EventClient::buffer_writecb(struct bufferevent *bev, void *arg)
{;}

//
// call back on error or connect
//
void EventClient::buffer_eventcb(struct bufferevent *bev, short what, void *arg)
{
  if( what & BEV_EVENT_ERROR )
    {
      //if( errno == ECONNREFUSED )
      //warn("client error %d %d", what, errno);
      //}
      
      // free event buffer
      (static_cast<EventClient*>(arg))->setConnected( false );
      (static_cast<EventClient*>(arg))->removeBufEvent();
    }
  else if(what & BEV_EVENT_CONNECTED)
    {
      //evutil_socket_t fd = bufferevent_getfd(bev);
      (static_cast<EventClient*>(arg))->setConnected( true );
    }
}


#endif // WITLIBEVENT2
