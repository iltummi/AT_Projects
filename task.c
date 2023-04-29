#include <stdio.h>
#include <pthread.h>
#include <time.h>


#include "task.h"
#include "data_structure.h"
#include "time_operations.h"


int task_create(void*(*task)(void *),int i, int period, int drel, int prio) {
	pthread_attr_t myatt;
	struct sched_param mypar;
	int tret;
	if (i > 3) return -1;
	tp[i].arg = i;
	tp[i].period = period;
	tp[i].deadline = drel;
	tp[i].priority = prio;
	tp[i].dmiss = 0;

	pthread_attr_init(&myatt);
	pthread_attr_setinheritsched(&myatt,PTHREAD_EXPLICIT_SCHED);
	pthread_attr_setschedpolicy(&myatt, SCHED_FIFO);
	mypar.sched_priority = tp[i].priority;
	pthread_attr_setschedparam(&myatt, &mypar);
	
	tret = pthread_create(&tid[i], &myatt, task,(void*)(&tp[i]));
	return tret;
}
int task_delete(int i){
	pthread_cancel(tid[i]);
	return 1;
}


int get_task_index(void* arg) {
	struct task_par *tpar;
	tpar = (struct task_par *)arg;
	return tpar->arg;
}

void set_activation(int i) {
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	time_copy(&(tp[i].at), t);
	time_copy(&(tp[i].dl), t);
	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].deadline);
}

int deadline_miss(int i) {
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	if (time_cmp(now, tp[i].dl) > 0) {
	tp[i].dmiss++;
	return 1;
	}
return 0;
}

void wait_for_activation(int i) {
	clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME, &(tp[i].at), NULL);
	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].period);
}



