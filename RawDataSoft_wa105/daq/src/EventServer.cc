//////////////////////////////////////////////////////////////////////////////////////////
//
//  Simple event dispatch server. It needs to be compiled with 
//  CFLAGS -lpthread -levent -levent_pthreads
//  The server spins a signle thread which does event_base_dispatch thus making the  
//  main thread non-blocking
//
//  Principal functions
//    - EventServer::Start() starts event server. The server blocking function runs
//                           in a separate thread
//    - EventServer::Stop() stops server thread and joins it with the calling thread
//    - EventServer::Send(const char* buf, size_t nb) buffer with number of bytes
//
//  Created: vgalymov@ipnl.inp2p3.fr 
//  Modified:
//
//////////////////////////////////////////////////////////////////////////////////////////

#ifdef WITHLIBEVENT2

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/socket.h>

#include "EventServer.h"
#include "LogMsg.h"

using namespace std;
using namespace dlardaq;

// set socket flags to non-blocking
int dlardaq::setnonblock(int &fd)
{
  int flags = fcntl(fd, F_GETFL);
  if( flags < 0 ) return flags;
  flags |= O_NONBLOCK;
  if (fcntl(fd, F_SETFL, flags) < 0) return -1;
  return 0;
}

// open a non-blocking listening socket
int dlardaq::openlistener(struct sockaddr_in &sin, int backlog)
{
  int fd;
  
  // create listener socket
  if( (fd = socket(AF_INET, SOCK_STREAM, 0)) < 0 )
    {
      msg_err<<"failed to create listening socket"<<endl;
      exit(1);
    }
  
  // bind it
  if( bind(fd, (struct sockaddr*)&sin, sizeof(sin)) < 0 )
    {
      msg_err<<"failed to bind listening socket on port "<<ntohs(sin.sin_port)<<endl;
      exit(1);
    }
  
  // mark it as passive socket
  if( listen(fd, MAXCLIENTS ) < 0)
    {
      msg_err<<"failed execute listen"<<endl;
      exit(1);
    }
  
  int opt = 1;
  if( setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt)) < 0 )
    {
      msg_err<<"failed to set socket options"<<endl;
      exit(1);
    }
  
  // set socket to be non-blocking
  if( setnonblock(fd) < 0 )
    {
      msg_err<<"failed to set socket to non-blocking"<<endl;
      exit(1);
    }

  return fd;
}


//
// ctor
//
EventServer::EventServer(uint16_t port, const char *ip)
{
  memset(&m_sin, 0, sizeof(m_sin));
  m_sin.sin_family = AF_INET;
  m_sin.sin_port   = htons(port);
  if(ip) inet_aton(ip, &m_sin.sin_addr);
  else m_sin.sin_addr.s_addr = htonl(INADDR_ANY);
  
  m_running = false;  
  Init();
}

//
// dtor
//
EventServer::~EventServer()
{
  if(m_running) Stop();
  
  // remove all clients
  while( !m_clients.empty() )
    RemoveClient(m_clients.begin()->second);
  
  event_free(ev_accept);
  event_base_free(ev_base);  
  
  close( m_listenfd );  
}

//
// Initializations
//
void EventServer::Init()
{
  m_listenfd = openlistener(m_sin, MAXCLIENTS);
  
  //
  evthread_use_pthreads();

  // create event base
  if( (ev_base = event_base_new()) == NULL )
    {
      close(m_listenfd);
      //err(1,"failed to create event base");
      msg_err<<"failed to create event base"<<endl;
      exit(1);
    }

  
  ev_accept = event_new( ev_base, m_listenfd, EV_READ | EV_PERSIST, 
			 EventServer::accept_cb, this );

  event_add(ev_accept, NULL);

  // handle interrupt signal events
  event_add(evsignal_new(ev_base, SIGPIPE, EventServer::signal_cb, this), NULL);
  //event_add(evsignal_new(ev_base, SIGHUP, EventServer::signal_cb, this), NULL);
  //event_add(evsignal_new(ev_base, SIGINT, EventServer::signal_cb, this), NULL);
}


//
// this function is non-blocking
//
void EventServer::Start()
{
  if(m_running) return;

  // start listener thread
  if( pthread_create(&m_run_thread, NULL, threadProcess, this) != 0 )
    {
      msg_err<<"Failed to create thread"<<endl;
      return;
    }
  
  m_running = true;
  
  //
  msg_info<<"Started data echo server"<<endl;
}


//
// this function should not block in principle
//
void EventServer::Stop()
{
  if(!m_running) return;
  m_running = false;
  stopThread();
  
  // join the dispatch thread
  pthread_join( m_run_thread, NULL);
  
  //
  msg_info<<"Stopped data echo server"<<endl;
}

//
//
//
void EventServer::runThread()
{
  //
  //event_add(ev_accept, NULL);
  event_base_dispatch(ev_base);
}

//
//
//
void EventServer::stopThread()
{
  //
  //event_del(ev_accept);
  event_base_loopexit(ev_base, NULL);  
}


//
// accept_conn_cb
//
void EventServer::accept_cb(evutil_socket_t fd, short what, void *arg)
{
  EventServer* obj = (EventServer*)arg;
  
  int cli_fd;
  struct sockaddr_in cli_addr;
  socklen_t cli_len = sizeof(cli_addr);

  cli_fd = accept(fd, (struct sockaddr *)&cli_addr, &cli_len);
  if (cli_fd < 0) 
    {
      msg_warn<<"Accept client failed"<<endl;
      return;
    }
  
  if(obj->GetClients() == MAXCLIENTS)
    {
      msg_warn<<"Cannot accept more clients"<<endl;
      close(cli_fd);
      return;
    }

  // set client socket to non-blocking mode
  if(setnonblock(cli_fd) < 0) 
    {
      msg_warn<<"Failed to set client socket to non-blocking"<<endl;
      close(cli_fd);
      return;
    }

  
  // create new client structure
  client_t* client = (client_t*)malloc(sizeof(*client));
  memset(client, 0, sizeof(*client));

  // client's fd
  client->fd = cli_fd;
  
  // add client
  obj->AddClient( client );
}


//
//
//
void EventServer::AddClient(client_t *client)
{
  std::map< int, client_t* >::iterator it;
  it = m_clients.find( client->fd );
  if( it!=m_clients.end() ) //should not happen
    {
      //warn("found client with the same fd");
      msg_warn<<"found client with the same fd"<<endl;
      RemoveClient( it->second );
    }

  // create bufferevent for this client socket
  client->base = ev_base;
  client->bev = bufferevent_socket_new( client->base, client->fd, 
					BEV_OPT_CLOSE_ON_FREE | BEV_OPT_THREADSAFE );
  bufferevent_setcb( client->bev, EventServer::buffer_readcb, EventServer::buffer_writecb, 
		     EventServer::buffer_eventcb, this );
  
  struct timeval sock_read_sec  = { RWTIMEOUT_SEC, 0 };
  struct timeval sock_write_sec = { RWTIMEOUT_SEC, 0 };
  bufferevent_set_timeouts( client->bev, &sock_read_sec, &sock_write_sec );
  //bufferevent_enable(bev, EV_WRITE | EV_READ); // WRITE is enabled by default

  // map it
  m_clients[client->fd] = client;
}


//
//
//
void EventServer::RemoveClient(client_t *client)
{
  RemoveClient( client->fd );
}

//
// remove hanging clients / close connections
//
void EventServer::RemoveClient(int fd)
{
  std::map< int, client_t* >::iterator it;
  it = m_clients.find( fd );
  if( it == m_clients.end() ) //should not happen
    {
      //warn("Could not find client in the map");
      msg_warn<<"Could not find client in the map"<<endl;
      return;
    }
  client_t* client = it->second;

  if(client->bev) bufferevent_free( client->bev );
  free(client);
  
  m_clients.erase( it );
}

//
// signal handling callback
//
void EventServer::signal_cb(evutil_socket_t sig, short what, void *arg)
{
  
  // signal handling of different interrupts
  
  switch( sig )
    {
    case SIGPIPE: 
      // ignore broken pipe and and let buffer_eventcb
      // take care of closing disconnected clients
      break;
      
    default:
      // warn("No default action on interrupt signal %d", (int)sig);
      msg_warn<<"No default action on interrupt signal "<<sig<<endl;
      break;
    }
}


//
// nothing to do for now
//
void EventServer::buffer_readcb(struct bufferevent *bev, void *arg)
{;}


//
// evbuffer_get_length(bufferevent_get_output(bev)) should be 0
// 
void EventServer::buffer_writecb(struct bufferevent *bev, void *arg)
{;}

//
// call back on error
//
void EventServer::buffer_eventcb(struct bufferevent *bev, short what, void *arg)
{
  //msg_warn<<"client error status "<<what<<endl;
  int fd = (int)bufferevent_getfd(bev);
  (static_cast<EventServer*>(arg))->RemoveClient( fd );
}


//
// send event buffer to connected clients
//
void EventServer::Send(const char *buf, size_t nb)
{
  msg_info<<"Send data "<<nb<<" bytes to "<<m_clients.size()<<" clients"<<endl;
  
  std::map< int, client_t* >::iterator it = m_clients.begin();
  for(;it!=m_clients.end();it++)
    {
      client_t* client = it->second;
      
      // write data to client's output buffer
      bufferevent_lock( client->bev );
      struct evbuffer *output = bufferevent_get_output(client->bev);
      if( evbuffer_get_length( output ) == 0)
	bufferevent_write(client->bev, buf, nb);
      bufferevent_unlock( client->bev );
    }
}


#endif // WITLIBEVENT2
