/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_SOCKET

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>

#ifdef _WIN32
  /* See http://stackoverflow.com/questions/12765743/getaddrinfo-on-win32 */
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501  /* Windows XP. */
  #endif
  #include <winsock2.h>
  #include <ws2tcpip.h>
#else
  /* Assume that any non-Windows platform uses POSIX-style sockets instead. */
  // #include <fcntl.h> // TODO: nothing appears to be used from fcntl.h?
  #include <unistd.h>
  #include <fcntl.h>
  #include <sys/types.h>
  #include <sys/socket.h>
  #include <netinet/in.h>
  #include <netdb.h>
  #include <sys/un.h>
  #include <arpa/inet.h>
  #include <errno.h>
  #include <sys/poll.h>
  #include <netinet/tcp.h>
#endif

#include "sock.h"


ssize_t recvn(int fd, char *vptr, size_t n, int flags) {
    size_t  nleft;
    ssize_t nread;
    char   *ptr;

    ptr = vptr;
    nleft = n;
    while (nleft > 0) {
	if ( (nread = recv(fd, ptr, nleft, flags)) < 0) {
	    if (errno == EINTR)
		nread = 0;	/* and call recv() again */
	    else
		return (-1);
	} else if (nread == 0)
	    break;		/* EOF */
	nleft -= nread;
	ptr += nread;
    }
    return (n - nleft); 	/* return >= 0 */
}

ssize_t sendn(int fd, const char* vptr, size_t n, int flags) {
    size_t nleft;
    ssize_t nwritten;
    const char* ptr;
    ptr = vptr;
    nleft = n;
    while(nleft > 0) {
        if ((nwritten = send(fd, ptr, nleft, flags)) <= 0) {
            if (nwritten < 0 && errno == EINTR)
                nwritten = 0;  /* and call send() again */
            else
                return -1;     /* error */
        }
        nleft -= nwritten;
        ptr += nwritten;
    }
    return n;
}

int
mbdyn_make_inet_socket(SOCKET* sock, struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int dobind, int *perrno)
{
	return mbdyn_make_inet_socket_type(sock, name, hostname, port, MBDYN_DEFAULT_SOCKET_TYPE, dobind, perrno);
}

int
mbdyn_host2inet_addr(struct sockaddr_in *name, const char *hostname, unsigned short int port, int socket_type, int *perrno)
{
	if (hostname) {
#if defined(HAVE_GETADDRINFO)
		char portbuf[sizeof("65535") + 1];
		struct addrinfo hints = { 0 }, *res = NULL;
		int rc;

		rc = snprintf(portbuf, sizeof(portbuf), "%d", (int)port);
		if (rc > STRLENOF("65535")) {
			return -4;
		}

		hints.ai_family = AF_INET;
		hints.ai_socktype = socket_type;
		rc = getaddrinfo(hostname, portbuf, &hints, &res);
		if (rc != 0) {
			*perrno = WSAGetLastError();
			return -3;
		}

		name->sin_addr = ((struct sockaddr_in *)res->ai_addr)->sin_addr;

		freeaddrinfo(res);

#elif defined(HAVE_GETHOSTBYNAME)
		struct hostent *hostinfo;

		/* TODO: use getnameinfo() if available */
		hostinfo = gethostbyname(hostname);
		if (hostinfo == NULL) {
			*perrno = h_errno;
			return -3;
		}

		name->sin_addr = *(struct in_addr *)hostinfo->h_addr;
#elif defined(HAVE_INET_ATON)
		struct in_addr addr;
		if (inet_aton(hostname, &addr) == 0) {
			*perrno = WSAGetLastError();
			return -3;
		}
		name->sin_addr = addr.s_addr;
#else
		return -3;
#endif
	} else {
		name->sin_addr.s_addr = htonl(INADDR_ANY);
	}

	return 0;
}

int
mbdyn_make_inet_socket_type(SOCKET* sock, struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int socket_type, int dobind, int *perrno)
{
	struct sockaddr_in tmpname = { 0 };

	if (name == NULL) {
		name = &tmpname;
	}

	if (perrno) {
		*perrno = 0;
	}

   	/* Give the socket a name. */
   	name->sin_family = AF_INET;
   	name->sin_port = htons(port);

	int rc = mbdyn_host2inet_addr(name, hostname, port, socket_type, perrno);
	if (rc != 0) {
		return rc;
	}

   	/* Create the socket. */
	*sock = socket(PF_INET, socket_type, 0);
   	//bool isblocking = true;
   	if (*sock == INVALID_SOCKET) {
		if (perrno) {
			*perrno = WSAGetLastError();
		}
		return -1;
   	}

   	/* disable Nagle's algorithm (only for TCP sockets) */
	int local_flag = 1;
	int result = 1;
	if (socket_type == SOCK_DGRAM) {
		result = setsockopt(*sock,           /* socket affected */
				IPPROTO_TCP,     /* set option at TCP level */
				TCP_NODELAY,     /* name of option */
				(char *) &local_flag,  /* the cast is historical cruft */
				sizeof(int));    /* length of option value */
	} 	

	if (result != 0) {
		fprintf(stderr, "Unable to disable Nagle's algorithm, sockets may be slow\n");
	}
	if (dobind) {
		rc = bind(*sock, (struct sockaddr *) name, sizeof(struct sockaddr_in));
		if (rc == SOCKET_ERROR) {
			if (perrno) {
				*perrno = WSAGetLastError();
			}
			return -2;
		}
   	}

   	return 0;
}



const char* sock_err_string (int err)
{
#ifdef _WIN32
    return winsock_err_string (err);
#else
    return strerror(err);
#endif /* _WIN32 */
}

#ifdef _WIN32
const char* winsock_err_string (int err)
{
    switch (err) {
    case 0:                  return "No error";
    case WSAEINTR:           return "Interrupted system call";
    case WSAEBADF:           return "Bad file number";
    case WSAEACCES:          return "Permission denied";
    case WSAEFAULT:          return "Bad address";
    case WSAEINVAL:          return "Invalid argument";
    case WSAEMFILE:          return "Too many open sockets";
    case WSAEWOULDBLOCK:     return "Operation would block";
    case WSAEINPROGRESS:     return "Operation now in progress";
    case WSAEALREADY:        return "Operation already in progress";
    case WSAENOTSOCK:        return "Socket operation on non-socket";
    case WSAEDESTADDRREQ:    return "Destination address required";
    case WSAEMSGSIZE:        return "Message too long";
    case WSAEPROTOTYPE:      return "Protocol wrong type for socket";
    case WSAENOPROTOOPT:     return "Bad protocol option";
    case WSAEPROTONOSUPPORT: return "Protocol not supported";
    case WSAESOCKTNOSUPPORT: return "Socket type not supported";
    case WSAEOPNOTSUPP:      return "Operation not supported on socket";
    case WSAEPFNOSUPPORT:    return "Protocol family not supported";
    case WSAEAFNOSUPPORT:    return "Address family not supported";
    case WSAEADDRINUSE:      return "Address already in use";
    case WSAEADDRNOTAVAIL:   return "Can't assign requested address";
    case WSAENETDOWN:        return "Network is down";
    case WSAENETUNREACH:     return "Network is unreachable";
    case WSAENETRESET:       return "Net connection reset";
    case WSAECONNABORTED:    return "Software caused connection abort";
    case WSAECONNRESET:      return "Connection reset by peer";
    case WSAENOBUFS:         return "No buffer space available";
    case WSAEISCONN:         return "Socket is already connected";
    case WSAENOTCONN:        return "Socket is not connected";
    case WSAESHUTDOWN:       return "Can't send after socket shutdown";
    case WSAETOOMANYREFS:    return "Too many references, can't splice";
    case WSAETIMEDOUT:       return "Connection timed out";
    case WSAECONNREFUSED:    return "Connection refused";
    case WSAELOOP:           return "Too many levels of symbolic links";
    case WSAENAMETOOLONG:    return "File name too long";
    case WSAEHOSTDOWN:       return "Host is down";
    case WSAEHOSTUNREACH:    return "No route to host";
    case WSAENOTEMPTY:       return "Directory not empty";
    case WSAEPROCLIM:        return "Too many processes";
    case WSAEUSERS:          return "Too many users";
    case WSAEDQUOT:          return "Disc quota exceeded";
    case WSAESTALE:          return "Stale NFS file handle";
    case WSAEREMOTE:         return "Too many levels of remote in path";
    case WSASYSNOTREADY:     return "Network system is unavailable";
    case WSAVERNOTSUPPORTED: return "Winsock version out of range";
    case WSANOTINITIALISED:  return "WSAStartup not yet called";
    case WSAEDISCON:         return "Graceful shutdown in progress";
    case WSAHOST_NOT_FOUND:  return "Host not found";
    case WSANO_DATA:         return "No host data of that type was found";
    default:                 return "Unknown winsock error";
    }
}

#else
/* errno cannot be used on windows for sockets, so we use
   WSAGetLastError but on non-windows define it to just return
   errno */
int WSAGetLastError(void)
{
    return errno;
}

int
mbdyn_make_named_socket(SOCKET* sock, struct sockaddr_un *name, const char *path,
	int dobind, int *perrno)
{
	return mbdyn_make_named_socket_type(sock, name, path, MBDYN_DEFAULT_SOCKET_TYPE, dobind, perrno);
}

int
mbdyn_make_named_socket_type(SOCKET *sock, struct sockaddr_un *name, const char *path,
	int socket_type, int dobind, int *perrno)
{
   	struct sockaddr_un tmpname = { 0 };
	socklen_t size;
	size_t pathlen;

   	*sock = INVALID_SOCKET;

	pathlen = strlen(path);
	if (pathlen >= sizeof(tmpname.sun_path)) {
		if (perrno) {
			*perrno = ENAMETOOLONG;
		}
		return -3;
	}

	if (name == NULL) {
		name = &tmpname;
	}

	if (perrno) {
		*perrno = 0;
	}

   	/* Create the socket. */
   	*sock = socket(PF_LOCAL, socket_type, 0);
   	if (*sock < 0) {
		if (perrno) {
			*perrno = errno;
		}
      		return -1;
   	}

   	/* Give the socket a name. */
   	name->sun_family = AF_LOCAL;
	// sizeof(name->sun_path)-1 because otherwise name->sun_path would not be null-terminated
	// for path longer than sizeof(name->sun_path)
   	strncpy(name->sun_path, path, sizeof(name->sun_path)-1);
#ifdef HAVE_OFFSETOF
	size = (offsetof(struct sockaddr_un, sun_path) + pathlen + 1);
#else /* HAVE_OFFSETOF */
	/* NOTE: not robust
	size = sizeof(struct sockaddr_un) + pathlen + 1 - sizeof(name->sun_path);
	*/
	/* perhaps this is better: */
	size = (void *)&tmpname.sun_path - (void *)&tmpname + pathlen + 1;
#endif /* !HAVE_OFFSETOF */

   	if (dobind) {
		int rc = bind(*sock, (struct sockaddr *)name, size);
		if (rc < 0) {
			if (perrno) {
				*perrno = errno;
			}

      			return -2;
		}
   	}

   	return 0;
}
#endif /* _WIN32 */

#endif /* USE_SOCKET */
