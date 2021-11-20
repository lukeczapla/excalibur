
#include "thread.h"


void *jthread_periodic_wrapper(void *jthread_ptr) {
	void *returnval = NULL;
	struct timeval tv;
	jthread *jt;
	
	jt = (jthread *)jthread_ptr;
	
	while (1) {
		pthread_mutex_lock(&jt->lock);
			if (jt->state != JTHREAD_RUN) {
				pthread_mutex_unlock(&jt->lock);
				break;
			}
			tv.tv_sec = jt->s_delay;
			tv.tv_usec = jt->mu_delay;
		pthread_mutex_unlock(&jt->lock);
		
		returnval = jt->func(jt);
		
		usleep((jt->s_delay*1000000) + jt->mu_delay);
	}
	
	pthread_mutex_lock(&jt->lock);
		jt->state = JTHREAD_STOPPED;
	pthread_mutex_unlock(&jt->lock);

	return returnval;
}


jthread *GetJThread(void *(*func)(void *), void *arg) {
	int result;
	jthread *jt;
	
	#ifdef DEBUG
		if (func == NULL) {
			printf("%s(): NULL func() passed\n", __func__);
			return NULL;
		}
	#endif
	
	jt = (jthread *)malloc(sizeof(jthread));
	if (jt == NULL) {
		printf("problem in %s(): malloc failed.\n", __func__);
		return NULL;
	}
	
	jt->func = func;
	jt->arg = arg;
	jt->s_delay = -1;
	jt->mu_delay = -1;
	jt->state = JTHREAD_STOPPED;

	result = pthread_mutex_init(&jt->lock, NULL);
	if (result != 0) {
		printf("problem in %s(): pthread_mutex_init: %s\n", __func__, strerror(result));
		free(jt);
		return NULL;
	}
	
	return jt;
}


void *StopJThread(jthread *jt, int wait_on_thread) {
	int result;
	void *returnval = NULL;
	
	#ifdef DEBUG
		if (jt == NULL) {
			printf("%s(): NULL jthread passed\n", __func__);
			return NULL;
		}
	#endif

	pthread_mutex_lock(&jt->lock);
		if (jt->state == JTHREAD_RUN)
			jt->state = JTHREAD_STOPPING;
		else {
			pthread_mutex_unlock(&jt->lock);
			return NULL;
		}
	pthread_mutex_unlock(&jt->lock);

	if (wait_on_thread == 1) {
		result = pthread_join(jt->thread, &returnval);
		if (result != 0) {
			printf("problem in %s(): pthread_cancel: %s\n", __func__, strerror(result));
		}
	}
	return returnval;
}


int RunJThread(jthread * jt) {
	int result;
	
	#ifdef DEBUG
		if (jt == NULL) {
			printf("%s(): NULL jthread passed\n", __func__);
			return -1;
		}
	#endif

	pthread_mutex_lock(&jt->lock);
		if(jt->state == JTHREAD_RUN)
		{
			pthread_mutex_unlock(&jt->lock);
			return 0;
		}
		jt->state = JTHREAD_RUN;
	pthread_mutex_unlock(&jt->lock);

	result = pthread_create(&jt->thread, NULL, jt->func, jt->arg);
	if (result != 0) return -1;
	return 0;
}


/*
	If you run this on a currently running periodic thread, it changes the periodic delay.
*/
int RunJThreadPeriodic(jthread * jt, int s_delay, int mu_delay) {
	int result;
	
	#ifdef DEBUG
		if (jt == NULL) {
			printf("%s(): NULL jthread passed\n", __func__);
			return -1;
		}
	#endif
	
	pthread_mutex_lock(&jt->lock);
		jt->s_delay = s_delay;
		jt->mu_delay = mu_delay;
		if (jt->state == JTHREAD_RUN) {
			pthread_mutex_unlock(&jt->lock);
			return 0;
		}
		jt->state = JTHREAD_RUN;
	pthread_mutex_unlock(&jt->lock);

	/*
		note; passed param is not jt->arg (as in RunJThread() ) because
		we also need some further info regarding the periodic delay which is
		stored in the jthread structure.
	*/
	result = pthread_create(&jt->thread, NULL, jthread_periodic_wrapper, jt);
	if (result != 0) return -1;
	return 0;
}


void *FreeJThread(jthread *jt) {
	void *returnval;
	
	#ifdef DEBUG
		if (jt == NULL) {
			printf("%s(): NULL jthread passed\n", __func__);
			return NULL;
		}
	#endif

	returnval = StopJThread(jt, 1);
	pthread_mutex_destroy(&jt->lock);
	free(jt);
	return returnval;
}

