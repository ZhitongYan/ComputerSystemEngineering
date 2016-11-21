// this is the lock server
// the lock client has a similar interface

#ifndef lock_server_h
#define lock_server_h

#include <string>
#include "lock_protocol.h"
#include "lock_client.h"
#include "rpc.h"

class lock_server {

 protected:
  int nacquire;

 private:
  enum lstate {FREE, OCCUPIED };
  pthread_mutex_t lock;
  pthread_cond_t cond;
  std::map<lock_protocol::lockid_t, lstate> lstat_map;
  std::map<lock_protocol::lockid_t, int> ltime_map;
  std::map<lock_protocol::lockid_t, int> lowner_map;

 public:
  lock_server();
  ~lock_server() {};
  lock_protocol::status stat(int clt, lock_protocol::lockid_t lid, int &);
  lock_protocol::status acquire(int clt, lock_protocol::lockid_t lid, int &);
  lock_protocol::status release(int clt, lock_protocol::lockid_t lid, int &);
};

#endif 







