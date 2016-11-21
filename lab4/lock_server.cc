// the lock server implementation

#include "lock_server.h"
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include <arpa/inet.h>

lock_server::lock_server():
  nacquire (0)
{
	pthread_mutex_init(&lock, NULL);
	pthread_cond_init(&cond, NULL);
	lstat_map.clear();
	ltime_map.clear();
	lowner_map.clear();	
}

lock_protocol::status
lock_server::stat(int clt, lock_protocol::lockid_t lid, int &r)
{
  lock_protocol::status ret = lock_protocol::OK;
  printf("stat request from clt %d\n", clt);
  pthread_mutex_lock(&lock);
  //r = nacquire;
  r = ltime_map[lid];
  pthread_mutex_unlock(&lock);
  return ret;
}

lock_protocol::status
lock_server::acquire(int clt, lock_protocol::lockid_t lid, int &r)
{
  lock_protocol::status ret = lock_protocol::OK;
	// Your lab4 code goes here
	printf("client %i wants to acquire %lld\n", clt, lid);
	if (clt < 0) {
		ret = lock_protocol::RPCERR;
		return ret;
	}
	pthread_mutex_lock(&lock);
	if (lstat_map.find(lid) != lstat_map.end()) {
		while (lstat_map[lid] == OCCUPIED) {
			printf("cond wait...\n");
			pthread_cond_wait(&cond, &lock);
		}
	}
	lstat_map[lid] = OCCUPIED;
	ltime_map[lid]++;
	lowner_map[lid] = clt;
	pthread_mutex_unlock(&lock);
  return ret;
}

lock_protocol::status
lock_server::release(int clt, lock_protocol::lockid_t lid, int &r)
{
  lock_protocol::status ret = lock_protocol::OK;
	// Your lab4 code goes here
	printf("client %i wants to release %lld\n", clt, lid);
	pthread_mutex_lock(&lock);
	if ((lstat_map.find(lid) == lstat_map.end()) || (lstat_map[lid] == FREE)) {
		ret = lock_protocol::NOENT;
		pthread_mutex_unlock(&lock);
		return ret;
	}
	if (lowner_map[lid] != clt) {
		ret = lock_protocol::RPCERR;
		pthread_mutex_unlock(&lock);
		return ret;
	}
	lstat_map[lid] = FREE;
	lowner_map[lid] = -1;
	pthread_cond_broadcast(&cond);
	pthread_mutex_unlock(&lock);
  return ret;
}
