#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <allegro.h>
#include <semaphore.h>

#include "data_structure.h"
#include "task.h"
#include "camera_unico.h"
	
void* scatta(void *arg){
		
		int ti;
		ti=get_task_index(arg);
		set_activation(ti);
		
		struct win finestra_locale;
		struct coord coord_palla_locale;
		
		coord_palla_locale.xnew = -1;
		coord_palla_locale.ynew = -1;
		
		int trovato;

		printf("sono scatta\n");
		
		while(1){
			
		system("python3 scatta.py");
		foto = load_bitmap("foto.bmp", NULL);
		draw_sprite(screen, foto, 0, 0);
		
		if (foto != NULL) { //foto trovata
			
			sem_wait(&sem_finestra);
			finestra_locale = finestra;
			sem_post(&sem_finestra);
			
			trovato = centroid(finestra_locale, &coord_palla_locale);
			
			if(trovato){
			sem_wait(&sem_coord);
			coord_palla = coord_palla_locale;
			printf("xnew = %d ynew = %d xold= %d yold = %d \n" , coord_palla.xnew, coord_palla.ynew, coord_palla.xold, coord_palla.yold);
			sem_post(&sem_coord);
			
			wait_for_activation(ti);
		    }
			else{
			printf("palla non trovata\n");
			sem_wait(&sem_coord);
			coord_palla.xnew = -1;
			coord_palla.ynew = -1;
			sem_post(&sem_coord);
			wait_for_activation(ti);}
			
			}
		else{
			printf("foto non trovata\n");
			wait_for_activation(ti); }
		}
}

int centroid(struct win w, struct coord *target){
		
		int c, x, y, xn, yn;
		long sum, xc, yc;
		sum=0;
		xc=0;
		yc=0;
		int r ,g , b;
		for (x=0; x<w.xsize; x++) {
			for (y=0; y<w.ysize; y++) {
				c = getpixel(foto, x+w.x, y+w.y);
				r = getr(c); g = getg(c); b = getb(c);
				if (r < 50 && g < 50 && b<50) { // search for target
					sum++;
					xc = xc + x;
					yc = yc + y;
				}
			}
		}
		if (sum > 0) { // target found
			xn = (xc + w.x)/sum;
			yn = (yc + w.y)/sum;
			aggiorna_posizione(target, xn, yn);
			return 1;
		}
			else return 0; // target not found
}

void aggiorna_posizione(struct coord *posizione,int x,int y){

	posizione->xold = posizione->xnew;
	posizione->yold = posizione->ynew;
	
	posizione->xnew = x;
	posizione->ynew = y;
}

