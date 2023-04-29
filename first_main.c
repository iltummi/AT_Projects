#include <stdio.h>
#include <sched.h>
#include <stdlib.h>
#include <allegro.h>
#include <semaphore.h>


#include "data_structure.h"
#include "task.h"
#include "prova_task.h"
#include "time_operations.h"
#include "motor.h"
#include "camera.h"
#include "predittore.h"


/*struct sched_param mypar;
struct task_par tp[NT];
pthread_t tid[NT];*/



int main(){
	
	lung=600;
	larg=490;
	T_camera=4000; 
	
	//motore
	inizializza_pin();
	//passi = 10;
	//dir = 1;
	//task_create(motor, 0, 1000, 1500, 90);
	/*
	//camera
	allegro_init();
	set_color_depth(16);
	set_gfx_mode(GFX_AUTODETECT_WINDOWED, 600, 490, 0, 0);
	
	sem_init(&sem_coord, 0, 1);
	sem_init(&sem_finestra, 0, 1);
	//sem_init(&sem_foto, 0, 1);
	sem_init(&sem_passi, 0, 1);
	sem_init(&sem_motor, 0, 0);
	
	finestra.x = 0;
	finestra.xsize = 600;
	finestra.y = 0;
	finestra.ysize = 490;
	
	task_create(scatta,0, T_camera, 2000 ,50) ;
	//foto = load_bitmap("foto.bmp", NULL);
	//task_create(trova_palla, 1, 2000, 1500, 20);
	task_create(predittore, 1, 4000, 1000, 20);
	*/
	while(1){

	printf("passi:");
	scanf("%d", &passi);	
	printf("dir:");
	scanf("%d", &dir);	
	task_create(motor, 0, 5, 150, 90);
	}


}
