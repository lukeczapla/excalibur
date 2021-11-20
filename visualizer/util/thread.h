#ifndef _VTHREAD_H
#define _VTHREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <pthread.h>

#define JTHREAD_NORMAL 		0
#define JTHREAD_PERIODIC	1

#define JTHREAD_RUN 		2
#define JTHREAD_STOPPING 	1
#define JTHREAD_STOPPED 	0


typedef struct {
	pthread_t thread;
	pthread_mutex_t lock;

	int type, state, s_delay, mu_delay;
	void *arg;
	void *(*func)(void *);
} jthread;


jthread *GetJThread(void *(*func)(void *), void * arg);
	
void *StopJThread(jthread *jt, int wait_on_thread);
	
int RunJThread(jthread *jt);
int RunJThreadPeriodic(jthread *jt, int s_delay, int mu_delay);
	
void *FreeJThread(jthread *jt);
	
#endif
