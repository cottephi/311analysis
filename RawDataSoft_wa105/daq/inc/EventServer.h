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
//    - EventServer::Send(const char* buf, size_t nb) send buffer of number of bytes
//
//  Created: vgalymov@ipnl.inp2p3.fr 
//  Modified: 
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __EVENTSERVER_H__
#define __EVENTSERVER_H__

#ifdef WITHLIBEVENT2 //only enabled by compiling with appropriate flags

#include <map>
#include <stdint.h>
#include <pthread.h>
#include <arpa/inet.h>
#include <event2/event.h>
#include <event2/thread.h>
#include <event2/buffer.h>
#include <event2/bufferevent.h>

#define MAXCLIENTS 2     // max number of clients this server accepts
#define RWTIMEOUT_SEC 2  // r/w timeout in seconds

//
namespace dlardaq
{
  int setnonblock(int &fd);
  int openlistener(struct sockaddr_in &sin, int backlog);  

  class EventServer
  {
  public:
    // can optionally bind listener socket on an alias IP 
    // default behaviour is to use INADDR_ANY for bind
    EventServer(uint16_t port, const char *ip = NULL);
    ~EventServer();
    
    // start the server
    void Start();

    // stop the server
    void Stop();

    // send data to clients
    void Send(const char* buf, size_t nb);
    
    // check if we are up
    bool IsRunning() const { return m_running; }

  private:
    // thread function
    static void* threadProcess(void *obj)
    {
      ((EventServer *)obj)->runThread();
      return NULL;
    }
    
    void Init();
    void runThread();
    void stopThread();
    
    // client object
    typedef struct client
    {
      int fd;                    // file descriptor for the client
      struct event_base*  base;  // should be this base
      struct bufferevent* bev;   // bufferevent for this client
    } client_t;
    

    size_t GetClients() const { return m_clients.size(); }
    void   AddClient(client_t* client); 
    void   RemoveClient(client_t* client);
    void   RemoveClient(int fd);
    
    // accept connection from clients
    static void accept_cb(evutil_socket_t fd, short what, void *arg);
    
    // callback for signals (sigpipe)
    static void signal_cb(evutil_socket_t sig, short what, void *arg); 
    
    // buffer callback functions
    static void buffer_readcb(struct bufferevent *bev, void *arg);
    static void buffer_writecb(struct bufferevent *bev, void *arg);
    static void buffer_eventcb(struct bufferevent *bev, short what, void *arg);

    // running 
    bool m_running;
    
    // number of clients
    std::map< int, client_t* > m_clients;

    // listener socket
    int m_listenfd;

    // server net stuff
    struct sockaddr_in m_sin;
    
    // run thread
    pthread_t m_run_thread;
        
    // lib event
    struct event*      ev_accept;
    struct event_base* ev_base;
  };  
};

#endif //WITLIBEVENT2

#endif //ifdef
