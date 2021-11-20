#ifndef _LJ_SOCKET_H
#define _LJ_SOCKET_H

#include <sys/socket.h>
#include <sys/select.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

typedef struct
{
	int socket;
	struct addrinfo *ai;
} jsocket;


jsocket *GetJSocket(const char *address, int port, int bind_socket);
void FreeJSocket(jsocket *js);
	
int ListenJSocket(jsocket *js, struct sockaddr_storage *s_in, socklen_t *s_len, int secs, int musecs);
	
int ReadJSocket(int socket, void *buf, int nbytes);
int ReadJSocketSized(int socket, void *buf, int maxbuf);

int WriteJSocket(int socket, void *buf, int nbytes);
int WriteJSocketSized(int socket, void *buf, int buflen);

int get_socketaddr_info(struct sockaddr_storage *ss, int *family, int *port, char *IP, int IP_max);
int GetJSocketInfo(jsocket *js, int *family, int *port, char *IP, int IP_max );


#endif

