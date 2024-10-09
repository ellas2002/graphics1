
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



int read_object(FILE *f, int obj_index){

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

int print_object (FILE *fout, int obj_index){

  fprintf(fout, "%d\n",numpoints[obj_index]) ;

  for (int i = 0 ; i < numpoints[obj_index] ; i++) {
    //12.6 ensures things are printed out with a
    //width of 12 charactesr includeing 6 decimal
    fprintf(fout, "%12.6lf %12.6lf\n",x[obj_index][i],y[obj_index][i]) ;
  }
  

  for (int i = 0 ; i < numpolys[obj_index] ; i++) {
    //integer shoud occupy at least three spaces
    fprintf(fout, "%3d    ",psize[obj_index][i]) ;

    for (int j = 0 ; j < psize[obj_index][i] ; j++) {
      //integer shoudl occupy at least two spaces
      fprintf(fout, "%2d ", con[obj_index][i][j]) ; }
    
    fprintf(fout, "\n") ;
  }


  for (int i = 0 ; i < numpolys[obj_index] ; i++) {
    fprintf(fout,"%lf %lf %lf\n",red[obj_index][i],grn[obj_index][i],blu[obj_index][i]) ;
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
  double t, temp, c, s;
  double xmin, xmax, ymin, ymax;
  double center_x, center_y;
  double dx, dy, temp_x, temp_y;

  xmin = xmax = x[objnum][0];
  ymin = ymax = y[objnum][0];
  

  for(int i = 0; i < numpoints[objnum]; i++){
    if(x[objnum][i] < xmin) xmin = x[objnum][i];
    if(x[objnum][i] > xmax) xmax = x[objnum][i];
    if(y[objnum][i] < ymin) ymin = y[objnum][i];
    if(y[objnum][i] > ymax) ymax = y[objnum][i];
  }
  
  //transiton center point to orgin (0,0)
  
  center_x = (xmin + xmax) / 2;
  center_y = (ymin + ymax) / 2;
  
  c = cos(radians);
  s = sin(radians);

  printf("Before rotation: First point (%.2f, %.2f)\n", x[objnum][0], y[objnum][0]);

  for (int i = 0; i < numpoints[objnum]; i++) {
       // Translate point to origin
    temp_x = x[objnum][i] - center_x;
    temp_y = y[objnum][i] - center_y;
        
        // Rotate point
    x[objnum][i] = temp_x * c - temp_y * s + center_x;
    y[objnum][i] = temp_x * s + temp_y * c + center_y;
    
  }
}




int main() {
    FILE *fin;
    int key, objnum = 0;
    char fname[100];
    double click_point[2];
    int num_objects = 0;
    double spin_angle = 15.0;

    while (num_objects < MAX_OBJECTS) {
        printf("Enter name of xy file (or 'q' to quit): ");
        scanf("%s", fname);

        if (fname[0] == 'q') break;

        fin = fopen(fname, "r");
        if (fin == NULL) {
            printf("Can't read file: %s\n", fname);
            continue;
        }

        read_object(fin, num_objects);
	center_object(num_objects);
	//spin_obj(num_objects, 2);
	
        //fclose(fin);
        num_objects++;
    }

    if (num_objects == 0) {
        printf("No objects loaded. Exiting.\n");
        return 0;
    }

    G_init_graphics(800, 800);

    // Draw the first object immediately
    G_clear();
    draw_object(objnum);

    while (true) {
        key = G_wait_key();
        printf("Key pressed: %c\n", key);

        if (key >= 'a' && key <= 'o') {
            int temp = key - 'a';
            if (temp < num_objects) {
                objnum = temp;
                G_rgb(0,0,0);
                G_clear();
                draw_object(objnum);
            }
        } 
        else if (key == 's' || key == 'S') {
            printf("Spinning object %d\n", objnum);
            spin_obj(objnum, spin_angle * M_PI / 180.0);
            G_rgb(0,0,0);
            G_clear();
            draw_object(objnum);
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



