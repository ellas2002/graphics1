

// open a sequence of .xy  files specified on the commmand line
// and draw them.

#include "FPToolkit.c"
#include <stdbool.h>
#include <math.h>
#define MAXPTS 59000
#define MAXPOLYS 57500
#define MAX_OBJECTS 100


int numpoints[MAX_OBJECTS] ;
int numpolys[MAX_OBJECTS] ;
double x[MAX_OBJECTS][MAXPTS] ;
double y[MAX_OBJECTS][MAXPTS] ;
int psize[MAX_OBJECTS][MAXPOLYS] ;
int con[MAX_OBJECTS][MAXPOLYS][20] ;
double red[MAX_OBJECTS][MAXPOLYS],grn[MAX_OBJECTS][MAXPOLYS],blu[MAX_OBJECTS][MAXPOLYS] ;



int read_object(int obj_index, FILE *f){

    // reads number of points from a file
    fscanf(f,"%d",&numpoints[obj_index]) ;

    //if the point info is more than maxpoints
    if (numpoints[obj_index] >= MAXPTS) {
      // need an extra for object centering
      printf("MAXPTS = %d :  exceeded.\n",MAXPTS) ;
      //kill
      exit(1) ;
    }

    for (int i = 0 ; i < numpoints[obj_index] ; i++) {
      //scans to see the x and y coordinates of the file
      fscanf(f,"%lf %lf",&x[obj_index][i],&y[obj_index][i]) ;
    }

    // connectivity info
    //reads number of polygons
    fscanf(f,"%d",&numpolys[obj_index]) ;
    //check to see if n of polygons is higher than MaxPoly
    if (numpolys[obj_index] > MAXPOLYS) {
      printf("MAXPOLYS = %d :  exceeded.\n",MAXPOLYS) ;
      //if so kill
      exit(1) ;
    }
    
    
    for (int i = 0 ; i < numpolys[obj_index] ; i++) {
      //scan for how many vertices the polygon has
      fscanf(f,"%d",&psize[obj_index][i]) ;
      for (int j = 0 ; j < psize[obj_index][i] ; j++) {
	//checking for the storage of indices of vertices
	//that make up each polygon in a 2d array
        fscanf(f,"%d",&con[obj_index][i][j]) ;
      } 
    } 

    
    // color info :
    for (int i = 0 ; i < numpolys[obj_index] ; i++) {
      fscanf(f,"%lf %lf %lf",&red[obj_index][i],&grn[obj_index][i],&blu[obj_index][i]) ;
    }    
}


int draw_object(int obj_index){
  double xp[100], yp[100] ;
  int np ;

  for (int i = 0 ; i < numpolys[obj_index] ; i++) {

    np = psize[obj_index][i] ; //get number of vertices in this polygon
    
    //get vertix coordinates for current polygon
    for (int j = 0 ; j < np ; j++) {
      int h = con[obj_index][i][j] ; // get index of vertex
      xp[j] = x[obj_index][h] ; //store the x coordinate
      yp[j] = y[obj_index][h] ; //store y coordiante
    }

    G_rgb(red[obj_index][i], grn[obj_index][i], blu[obj_index][i]) ;
    
    G_fill_polygon(xp, yp, np) ;
  }
}


void center_object(int objnum){
  //find center of bounding box string min x, min y, max x, maxy
  double xmin, xmax, ymin, ymax;
  double center_x, center_y;
  double dx, dy;

  xmin = xmax = x[objnum][0];
  ymin = ymax = y[objnum][0];
  
  //bounding box
  for(int i = 0; i < numpoints[objnum]; i++){
    if(x[objnum][i] < xmin) xmin = x[objnum][i];
    if(x[objnum][i] > xmax) xmax = x[objnum][i];
    if(y[objnum][i] < ymin) ymin = y[objnum][i];
    if(y[objnum][i] > ymax) ymax = y[objnum][i];
  }
  
  //transiton center point to orgin (0,0)
  center_x = (xmin + xmax) / 2;
  center_y = (ymin + ymax) / 2;
  
  //then to center screen (400,400)
  dx = 400 - center_x;
  dy = 400 - center_y;
  

  //apply some translate max to all x[], and y[] in the shape to move whole thing
  for (int i = 0; i < numpoints[objnum]; i++) {
        x[objnum][i] += dx;
        y[objnum][i] += dy;
    }
}

	    


void spin_obj(int objnum, double radians){
  double t = radians *M_PI/180;

  printf("Before rotation: First point (%.2f, %.2f)\n", x[objnum][0], y[objnum][0]);

  for (int i = 0; i < numpoints[objnum]; i++) {
    double temp = cos(t) * (x[objnum][i]- 400) - sin(t) * (y[objnum][i] - 400) + 400;
    y[objnum][i] = (x[objnum][i]- 400) * sin(t) + cos(t) *(y[objnum][i]- 400) + 400;
    x[objnum][i] = temp;
  }
}


int main(int argc, char ** argv) {
    FILE *fin;
    int key, obj_num;
    char fname[100];
    int num_objects = 0;
    double spin_angle = 15.0;

    for(int i = 0; i < argc; i++){
      if(i + 1 < argc){
	printf("argv[%d] is %s\n", i + 1, argv[i + 1]);
      }
    }

    num_objects = argc -  1;

    for(obj_num = 1; obj_num < num_objects + 1; obj_num++){
       fin = fopen(argv[obj_num], "r") ;

    
       if (fin == NULL) {
	 printf("can't read file, %s\n", argv[obj_num]) ;
	 exit(1) ;
       }
       
       center_object(obj_num - 1);
       read_object(obj_num - 1, fin);
       fclose(fin);
    }

    G_init_graphics(800, 800);
    G_clear();


    while (true) {
        key = G_wait_key();
        printf("Key pressed: %c\n", key);

        if (key >= 'a' && key <= 'o') {
            int temp = key - 'a';
            if (temp < num_objects) {
                obj_num = temp;
		G_rgb(0,0,0);
                G_clear();
		center_object(obj_num);
                draw_object(obj_num);
            }
        } 
        else if (key == 's' || key == 'S') {
            printf("Spinning object %d\n", obj_num);
            spin_obj(obj_num, spin_angle);
            G_rgb(0,0,0);
            G_clear();
            draw_object(obj_num);
            printf("Spin complete\n");
        }
        else if (key == 'q' || key == 'Q') {
            printf("Quitting...\n");
            break;
        }
    }

    G_close();
    return 0;
}



