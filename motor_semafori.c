#include <wiringPi.h>
#include <semaphore.h>
#include "task.h"
#include "motor.h"
#include "data_structure.h"

#define ENABLE 0
#define DIR 1
#define STEP 3

int up;

void* step(int passi){  //genera i passi del motore
	int up1;
	int ti;
	ti=get_task_index(step);
	set_activation(ti);
	
	while(1){
		if(up==1){
			digitalWrite(STEP, HIGH);
			up1=0;
			contastep++;
		}
		else if (up==0){
			digitalWrite(STEP, LOW);
			up1=1;
			contastep++;
			}
		up=up1;
		if(contastep==passi)
			sem_post(&sem_step);
			
		wait_for_activation(ti);
		}
	}



int motor(int passi, int dir){
	
	wiringPiSetup();
	pinMode(ENABLE, OUTPUT);
	pinMode(STEP, OUTPUT);
	pinMode(DIR, OUTPUT);
	
	digitalWrite(ENABLE, HIGH);
	
	if(dir==1)
		digitalWrite(STEP, HIGH);
	else 
		digitalWrite(STEP, LOW);
		
	int up=1;
	int contastep=0;

	
	
	task_create(step(passi) ,1, 1000000, 1000 , 20);
	//aspetta che siano fatti tutti i passi
	sem_wait(&sem_step);
	
	task_delete(1);
	digitalWrite(ENABLE, LOW);
	return 1;
}
	
	

	











