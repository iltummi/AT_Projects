
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <allegro.h>
#include <semaphore.h>

#include "data_structure.h"
#include "task.h"
	
void* scatta(void *arg){
		
		int ti;
		ti=get_task_index(arg);
		set_activation(ti);
		
		BITMAP *foto_locale;
		printf("sono scatta\n");
	    while(1){
			
		system("python3 scatta.py");
		printf("sorridi\n");
		
		foto_locale = load_bitmap("foto.bmp", NULL);
		
		if (foto_locale != NULL) { //foto trovata
			sem_wait(&sem_foto);
			foto = foto_locale;
			//blit(foto_locale, foto, 0, 0, 0, 0, foto_locale->w, foto_locale->h);
			printf("foto trovata e caricata\n");
			//draw_sprite(screen, foto, 0, 0);
			sem_post(&sem_foto);
			
			wait_for_activation(ti); 
			}
		else{
			printf("foto non trovata\n");
			wait_for_activation(ti); }
		}
}

int centroid(struct win w, struct coord *target){
		BITMAP *foto_locale2;
		//foto_locale = create_bitmap(640, 480);
		
		sem_wait(&sem_foto);
		foto_locale2 = foto;
		//blit(foto, screen, 0, 0, 0, 0, foto->w, foto->h);
		
		//printf("prima\n");
		//draw_sprite(screen, foto, 0, 0);
		//printf("dopo\n");
		sem_post(&sem_foto);
		//printf("dopo i semafori\n");
		//draw_sprite(screen, foto_locale, 0, 0);	
		//printf("dopo screen\n");
		//blit(foto_locale, screen, 0, 0, 0, 0, foto_locale->w, foto_locale->h);
		
		int c, x, y;
		long sum, xc, yc;
		sum=0;
		xc=0;
		yc=0;

		for (x=0; x<w.xsize; x++) {
			for (y=0; y<w.ysize; y++) {
				printf("prima\n");
				c = getpixel(foto_locale2, x+w.x, y+w.y);
				printf("dopo\n");
				if (c == 0) { // search for target
					sum++;
					xc = xc + x;
					yc = yc + y;
				}
			}
		}
		destroy_bitmap(foto_locale2);
		if (sum > 0) { // target found
			target->x = (xc + w.x)/sum;
			target->y = (yc + w.y)/sum;
			return 1;
		}
			else return 0; // target not found
}

void* trova_palla(void *arg){
		
		int ti;
		ti=get_task_index(arg);
		set_activation(ti);
		
		printf("sono trova_palla\n");
	    
	    struct win finestra_locale;
	    struct coord coord_palla_locale;
	    int trovato;
	    
	    while(1){
			sem_wait(&sem_finestra);
			finestra_locale = finestra;
			sem_post(&sem_finestra);
			
			trovato = centroid(finestra_locale, &coord_palla_locale);
			
			if(trovato){
			sem_wait(&sem_coord);
			coord_palla = coord_palla_locale;
			printf("x = %d y = %d \n", coord_palla.x, coord_palla.y);
			sem_post(&sem_coord);
			
			wait_for_activation(ti);
		    }
			else{
			printf("palla non trovata\n");
			wait_for_activation(ti);}
		}
}
