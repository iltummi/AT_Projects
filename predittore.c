#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <allegro.h>
#include <semaphore.h>
#include <string.h>

#include "data_structure.h"
#include "task.h"
#include "motor.h"
#include "predittore.h"

float abs_float(float numero){
if (numero < 0){
numero = -numero;
return numero;}
else 
return numero;
}

void vai_a(struct coord posizione){
		float vx, vy; //velocità lungo i 2 assi
		float tempo;		
		int goalk;
		float goalk_float;
		printf("p.xnew = %d p.ynew = %d p.xold= %d p.yold = %d \n" , posizione.xnew, posizione.ynew, posizione.xold, posizione.yold);
				
		vx =(float) (posizione.xnew - posizione.xold)/T_camera;
		vy = (float)  (posizione.ynew- posizione.yold)/T_camera;
		vx = abs_float(vx);
		tempo = (float) (lung - posizione.xnew)/vx;
		printf("vx = %2.3f, vy = %2.3f, tempo = %2.3f \n", vx, vy, tempo);
		goalk_float = posizione.ynew + (vy * tempo);
		goalk= (int) goalk_float;
		
		pos.old = pos.new;
		pos.new = goalk;	
		printf("pos.old= %d \n pos.new= %d \n", pos.old, pos.new);	
		}
		
void muoviti_di(){
		int distanza;
		float passi_float;
		distanza = pos.new - pos.old;
		
			
		passi_float = distanza/0.833; //0.833 per ruota con 120 denti
		sem_wait(&sem_passi);	
		passi = (int) passi_float;
		printf("questi sono i passi: %d \n", passi);
		sem_post(&sem_passi);
	}
	
void* predittore(void *arg){
	
	int ti;
		ti=get_task_index(arg);
		set_activation(ti);
	
	int dd;
	struct coord posizione;
	 
	 printf("sono predittore\n");
	while(1){
		
		sem_wait(&sem_coord);
		posizione = coord_palla;
		sem_post(&sem_coord);
		
		if(posizione.xnew != -1 || posizione.ynew != -1 ||posizione.xold != -1 || posizione.yold != -1){
			
			vai_a(posizione);
			
			muoviti_di();
			
			dd = task_create(motor, 2, 1000, 1000, 90);
			printf(" dd = %d\n" , dd);
			
			sem_wait(&sem_motor);
			wait_for_activation(ti);		
			}
		
		else 
			wait_for_activation(ti);
		}
}

			
