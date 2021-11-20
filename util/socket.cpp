#include "socket.h"

/*
	Prints an error message, closes the socket and frees the memory before returning NULL.
*/
#define JSOCK_ERROR(sock,call) \
{ \
	FreeJSocket( sock ); \
	return NULL; \
}



/*
	Gets basic info from sockaddr_storage structure.
*/
int get_socketaddr_info(struct sockaddr_storage *ss, int *family, int *port, char *IP, int IP_max) {

	struct sockaddr_in *sin4;
	struct sockaddr_in6 *sin6;
	
#ifdef DEBUG
	if (ss == NULL) return 0;
	if (family == NULL) return 0;
	if (port == NULL) return 0;
	if (IP == NULL) return 0;
	if (IP_max < 1) return 0;
#endif
	
	*family = ss->ss_family;
	
	if (ss->ss_family == AF_INET) {
		sin4 = (struct sockaddr_in *)ss;
		*port = sin4->sin_port;
		inet_ntop(sin4->sin_family, &sin4->sin_addr, IP, IP_max);
	}
	else {
		sin6 = (struct sockaddr_in6 *) ss;
		*port = sin6->sin6_port;
		inet_ntop(sin6->sin6_family, &sin6->sin6_addr, IP, IP_max);
	}
	
	return 1;
}


/*
	Just a wrapper for the above.
*/
int GetJSocketInfo(jsocket *js, int *family, int *port, char *IP, int IP_max) {
#ifdef DEBUG
	if( js == NULL ) return 0;
	if( family == NULL ) return 0;
	if( port == NULL ) return 0;
	if( IP == NULL ) return 0;
	if( IP_max < 1 ) return 0;
#endif

	return get_socketaddr_info( (struct sockaddr_storage *) js->ai->ai_addr, family, port, IP, IP_max );
}



/*
	In:	address is a string containing either an ip address in a.b.c.d form or a host name. If a NULL pointer passed, defaults to "localhost".
		port is the port to bind or connect on (will be converted into network byte order automatically)
		bind_socket is 0 where connect() should be called, else calls bind().
		
	Out:	valid jsocket *, or NULL.
*/
jsocket *GetJSocket(const char *address, int port, int bind_socket) {

	int i;
	struct addrinfo hints;
	char buffer[1024];
	jsocket * new_jsocket;
	
//	if (address == NULL) address = "localhost"; // don't need for getaddrinfo, as NULL = this host by default.

	/*
		Get socket structure memory, before we do anything.
	*/
	new_jsocket = (jsocket *)malloc( sizeof(jsocket) );
	if( new_jsocket == NULL )
	{
		printf( "%s(): unable to allocate memory for jsocket structure\n", __func__ );
		return NULL;
	}
		
	/*
		Set up some hints, specifically that we want IPv4 datagrams (IPv6 doesn't support broadcast datagrams).
	*/
	memset( &hints, 0, sizeof(hints) );
	hints.ai_family = AF_INET;
	/*
		For some reason, AF_UNSPEC won't work with numeric IPv6 address passed as "address";
		the connection always fails, with "No route to host" given.
	*/
//	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	/*
		Get the appropriate info regarding the destination of the broadcast.
	*/
	sprintf( buffer, "%d", port );
	i = getaddrinfo( address, buffer, &hints, &new_jsocket->ai );
	if( i != 0 ) JSOCK_ERROR( new_jsocket, gai_strerror(i) );


	/*
		Get a socket to send the info through
	*/
	new_jsocket->socket = socket( new_jsocket->ai->ai_family, new_jsocket->ai->ai_socktype, new_jsocket->ai->ai_protocol );
	if( new_jsocket->socket == -1 ) JSOCK_ERROR( new_jsocket, "socket()" );
	
	/*
		Ensure address is reusable immediately, so as not to tie up any resources on exit etc.
	*/
	i = 1;
	if( setsockopt( new_jsocket->socket, SOL_SOCKET, SO_REUSEADDR, &i, sizeof(i) ) == -1 ) JSOCK_ERROR( new_jsocket, "setsockopt()" );


	/*
		If we should connect(), attempt to connect.
		Else, attempt to bind() and listen().
	*/
	if( bind_socket == 0 )
	{
		if( connect( new_jsocket->socket, (struct sockaddr *) new_jsocket->ai->ai_addr, new_jsocket->ai->ai_addrlen ) == -1 ) JSOCK_ERROR( new_jsocket, "connect()" );
	}
	else
	{
		if( bind( new_jsocket->socket, (struct sockaddr *) new_jsocket->ai->ai_addr, new_jsocket->ai->ai_addrlen ) == -1 ) JSOCK_ERROR( new_jsocket, "bind()" );
		if( listen( new_jsocket->socket, 1 ) == -1 ) JSOCK_ERROR( new_jsocket, "listen()" );
	}

	/*
		Print a bit of info, for debug purposes.
	*/
	#ifdef DEBUG
			int j;
			GetJSocketInfo( new_jsocket, &i, &j, buffer, 1023 );
		printf( "Socket %s:%d (via %d) %s\n", buffer, port, j, bind_socket == 0 ? "connected" : "bound" );
	#endif

	return new_jsocket;
}
/*
	Simple - closes the socket and frees the memory.
*/
void FreeJSocket( jsocket * js )
{
	#ifdef DEBUG
		if( js == NULL )
		{
			// don't call JSOCK_ERROR, as calls this! recursive apocalypse!
			printf( "NULL jsockst pointer in FreeJSocket()\n" );
			return;
		}
	#endif

	close( js->socket );
	freeaddrinfo( js->ai );
	free( js );
}
/*
	In:	valid socket structure, sockaddr_in and socklen_t pointers
		secs and musecs are the seconds and microseconds to use in the nonblocking select() timeout. Where either is -1, blocking select() used.
	
	Out:	returns -1 on error, 0 on timeout, or otherwise a socket file descriptor.
			s_in and s_len will be filled with data from accept() call, if return > 0.
			
	Should really check for NULL pointers etc.
*/
int ListenJSocket( jsocket * js, struct sockaddr_storage * ss, socklen_t * s_len, int secs, int musecs )
{
	int rv;
	struct timeval tv;
	fd_set sockets;
	
	#ifdef DEBUG
		if( js == NULL ) JSOCK_ERROR( js, "NULL jsocket passed" );
		if( ss == NULL ) JSOCK_ERROR( js, "NULL sockaddr passed" );
		if( s_len == NULL ) JSOCK_ERROR( js, "NULL socklen_t passed" );
	#endif

	tv.tv_sec = secs;
	tv.tv_usec = musecs;
	
	FD_ZERO( &sockets );
	FD_SET( js->socket, &sockets );
	
	if( secs < 0 || musecs < 0 ) rv = select( js->socket+1, &sockets, NULL, NULL, NULL ); // blocking
	else rv = select( js->socket+1, &sockets, NULL, NULL, &tv ); // nonblocking, with timeout.
	
	if( rv > 0 && FD_ISSET( js->socket, &sockets ) )
	{
		*s_len = sizeof( struct sockaddr_storage );
		rv = accept( js->socket, (struct sockaddr *) ss, s_len );
		if( rv == -1 ) printf( "%s(): error %d in accept(): %s\n", __func__, errno, strerror(errno) );
	}
	if( rv == -1 ) printf( "%s(): error %d in select(): %s\n", __func__, errno, strerror(errno) );
	
	return rv;
}


/*
	Read and write don't actually need a jsocket structure, they work on the standard socket file descriptors.
	This is handy as we can use them on standard socket file descriptors, such as those from accept() etc,
	as the methods wrap the recv() and send() calls in error testing.
*/
int ReadJSocket(int socket, void *buf, int maxbuf) {
	int bytes_read, total_bytes_read, bytes_remaining;
	int32_t len;
	char *cptr;

#ifdef DEBUG
		if( buf == NULL )
		{
			printf( "%s(): NULL buffer passed, ignoring.\n", __func__ );
			return -1;
		}
		if( maxbuf < 1 )
		{
			printf( "%s(): maxbuf (%d) < 1, ignoring.\n", __func__, maxbuf );
			return -1;
		}
#endif

	/*
		Read data size
	*/
	cptr = (char *) &len;
	total_bytes_read = 0;
	bytes_remaining = sizeof(len);

	while( bytes_remaining > 0 )
	{
		bytes_read = recv( socket, &cptr[total_bytes_read], bytes_remaining, 0 );
		if( bytes_read == -1 ) { printf( "%s(): error reading data length via recv(): %s\n", __func__, strerror(errno) ); return -1; }
		if( bytes_read == 0 )  { printf( "%s(): error in reading data length via recv(), client closed socket\n", __func__ ); return -1; }

		bytes_remaining -= bytes_read;
		total_bytes_read += bytes_read;
	}
	
	len = ntohl(len);
	if (maxbuf < len) {
		printf("%s(): maxbuf (%d) < len (%d), abandoning.\n", __func__, maxbuf, len );
		return -1;
	}
	
	/*
		Read data
	*/
	cptr = (char *)buf;
	total_bytes_read = 0;
	bytes_remaining = len;
	
	while (bytes_remaining > 0) {
		bytes_read = recv(socket, &cptr[total_bytes_read], bytes_remaining, 0);
		if (bytes_read == -1) { printf( "%s(): error in recv(): %s\n", __func__, strerror(errno) ); return -1; }
		if (bytes_read == 0)  { printf( "%s(): error in recv(), client closed socket\n", __func__ ); return -1; }

		bytes_remaining -= bytes_read;
		total_bytes_read += bytes_read;
	}
	
	return total_bytes_read;
}



int WriteJSocket (int socket, void *buf, int buflen) {
	int bytes_written, total_bytes_written, bytes_remaining;
	int32_t len;
	char * cptr;

	#ifdef DEBUG
		if( buf == NULL )
		{
			printf( "%s(): NULL buffer passed, ignoring.\n", __func__ );
			return -1;
		}
		if( buflen < 1 )
		{
			printf( "%s(): buflen (%d) < 1, ignoring.\n", __func__, buflen );
			return -1;
		}
	#endif
	
	len = htonl( buflen );

	/*
		Write data size
	*/
	cptr = (char *) &len;
	total_bytes_written = 0;
	bytes_remaining = sizeof(len);

	while( bytes_remaining > 0 )
	{
		bytes_written = send( socket, &cptr[total_bytes_written], bytes_remaining, 0 );
		if( bytes_written == -1 ) { printf( "%s(): error writing data length via send(): %s\n", __func__, strerror(errno) ); return -1; }
		// send() doesn't return 0 for client closed socket; it's always -1 for error
		if( bytes_written < (int)sizeof(int32_t) ) { printf( "%s(): error in writing data length via send(), couldn't write %lu bytes for data length - got %d\n", __func__, sizeof(int32_t), bytes_written ); return -1; }

		bytes_remaining -= bytes_written;
		total_bytes_written += bytes_written;
	}
	
	/*
		Write data
	*/
	cptr = (char *)buf;
	total_bytes_written = 0;
	bytes_remaining = buflen;
	
	while( bytes_remaining > 0 )
	{
		bytes_written = send( socket, &cptr[total_bytes_written], bytes_remaining, 0 );
		if( bytes_written == -1 ) { printf( "%s(): error in send(): %s\n", __func__, strerror(errno) ); return -1; }

		bytes_remaining -= bytes_written;
		total_bytes_written += bytes_written;
	}

	return total_bytes_written;
}
