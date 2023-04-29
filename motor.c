#include <wiringPi.h>
#include <semaphore.h>


#include "task.h"
#include "motor.h"
#include "data_structure.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void inizializza_pin(){
	wiringPiSetup();
	pinMode(ENABLE, OUTPUT);
	pinMode(STEP, OUTPUT);
	pinMode(DIR, OUTPUT);
}

void* motor(void *arg){
	
	int ti;
	ti=get_task_index(arg);
	set_activation(ti);

	printf("hello motor \n");
	
	digitalWrite(ENABLE, LOW);//attiva motore

	int n_passi;
	
	//sem_wait(&sem_passi);
	n_passi= 2*passi; //un passo per ogni passaggio da 0 a 1
	//sem_post(&sem_passi);
	
	if (n_passi > 0)
			digitalWrite(DIR, HIGH);
	else
			digitalWrite(DIR, LOW);

	int up=1;
	int up1;

	int i=0;
	n_passi = abs(n_passi);
	
	for(i; i<n_passi; i++){
		if(up == 1){
			digitalWrite(STEP, HIGH);
			up1=0;
			printf("1\n");
		}
		else{
			digitalWrite(STEP, LOW);
			up1=1;
			printf("0\n");
			}
		up=up1;
		wait_for_activation(ti);
		}
	printf("finito \n");	
	digitalWrite(STEP, LOW);
	digitalWrite(ENABLE, LOW);
	//sem_post(&sem_motor);
}
	
	

	











