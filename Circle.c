 #include "FPToolkit.c"

int main()
{
  double swidth, sheight, x, y, i, t;
	swidth = 800; sheight = 800;
	G_init_graphics(swidth,sheight);

	G_rgb(0.3,0.3,0.3);
	G_clear();


	G_rgb(1.0,1.0, 0.0);
	G_circle(400,400,300);

	
	i = 0;
	while(i < 16){
	  t = i * 2 *M_PI/16;
	  x = cos(t)* 300 + 400;
	  y = sin(t)* 300 + 400;

	  G_line(400,400,x,y);
	  
	  i++;

	}

	int key;
	key = G_wait_key();

}
