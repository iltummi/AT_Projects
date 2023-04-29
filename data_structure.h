#include <time.h>
#include <stdio.h>
#include <pthread.h>
#include <allegro.h>
#include <semaphore.h>

#define NT 3 // n task

#define ENABLE 3
#define DIR 1
#define STEP 0

//Misure Tavola
int lung;
int larg;

int passi;
int dir;

BITMAP *foto; // pointer to bitmap
//foto = create_bitmap(640, 480);

//semafori camera.c
//sem_t sem_foto;
sem_t sem_coord;
sem_t sem_finestra;
sem_t sem_passi;
sem_t sem_motor;

//periodo task camera
int T_camera; 

extern struct posiz{
int old;
int new; 
	};
	
struct posiz pos;

extern struct win {
int x;
int y;
int xsize;
int ysize;
};

extern struct coord {
int xold;
int yold;
int xnew;
int ynew;
};

extern struct task_par{
int arg; /* task argument */
long wcet; /* in microseconds */
int period; /* in milliseconds */
int deadline; /* relative (ms) */
int priority; /* in [0,99] */
int dmiss; /* no. of misses */
struct timespec at; /* next activ. time */
struct timespec dl; /* abs. deadline */ };


struct sched_param	mypar;
struct task_par		tp[NT];
pthread_t			tid[NT];

struct win finestra;
struct coord coord_palla;
