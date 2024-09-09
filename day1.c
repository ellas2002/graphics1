#include "FPToolkit.c"

int main()
{
	int swidth, sheight;

	swidth = 800;
	sheight = 800;
	G_init_graphics (swidth, sheight);

	G_rgb (0.3,0.3,0.3);
	G_clear ();


	G_rgb(1,1,1);
	G_circle(400,400,300);


	int rot = 0;
	while (1){	
		G_rgb (0.3,0.3,0.3);
		G_clear ();
		G_rgb(1,1,1);
		G_circle(400,400,300);
		int i = 0;
			while (i<16){
				G_line(400,400,400+300*(cos(i*2*M_PI/16+rot)),400+300*(sin(i*2*M_PI/16+rot)));
				i++;
			}
		rot += 0.5;
		int key;
		key = G_wait_key();
	}
	int key;
	key = G_wait_key();
}
