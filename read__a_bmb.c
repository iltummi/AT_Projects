// from the terminal
// gcc read_a_bmb.c -o read_a_bmb `allegro-config --libs`
// ./read_a_bmb

#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <allegro.h>

int main(){
	
	allegro_init();
	set_color_depth(8);
	set_gfx_mode(GFX_AUTODETECT_WINDOWED, 400, 300, 0, 0);
	install_keyboard();
	
	
	BITMAP *fish; // pointer to bitmap
	int r, g, b, color;
	
	int i=0;
	
	for(i; i<30; i++){
		
		system("python3 bmb_try.py");
		fish = load_bitmap("snapshot.bmp", NULL);
		if (fish == NULL) {
			printf("file not found\n");
			exit(1);
		}
		color = getpixel(fish, 50, 250);
		r = getr(color);
		g = getg(color);
		b = getb(color);
		printf("color is %d, %d, %d, %d\n", color, r, g, b);
	
	}
	
	readkey();
	return 0;
}
