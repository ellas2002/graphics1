#include "FPToolkit.c"
#include "M2d_matrix_tools.c"
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
double matrix[3][3] = {0};


//reads the file in order to extract info and output the file
int read_object(int obj_index, FILE *f){

  //scan in number of points in object
  fscanf(f,"%d",&numpoints[obj_index]) ;

  //if the number of points is more than the alloted amount: kill
  if (numpoints[obj_index] >= MAXPTS) {
      printf("MAXPTS = %d :  exceeded.\n",MAXPTS) ;
      exit(1) ;
  }

  //scan in the x and y points
  for (int i = 0 ; i < numpoints[obj_index] ; i++) {
      fscanf(f,"%lf %lf",&x[obj_index][i],&y[obj_index][i]) ;
  }

  //scan in number of polygons
  fscanf(f,"%d",&numpolys[obj_index]) ;

  //if the number execeds alloted amount: kill
  if (numpolys[obj_index] > MAXPOLYS) {
      printf("MAXPOLYS = %d :  exceeded.\n",MAXPOLYS) ;
      exit(1) ;
   }
    
  //scan in polygon vertices (psize) and the storage of the indices of vertices (con)
  for (int i = 0 ; i < numpolys[obj_index] ; i++) {
     fscanf(f,"%d",&psize[obj_index][i]) ;
     for (int j = 0 ; j < psize[obj_index][i] ; j++) {
        fscanf(f,"%d",&con[obj_index][i][j]) ;
      } 
   } 
    
    // color info :
   for (int i = 0 ; i < numpolys[obj_index] ; i++) {
      fscanf(f,"%lf %lf %lf",&red[obj_index][i],&grn[obj_index][i],&blu[obj_index][i]) ;
   }    
}


//draws the object in the scanned file
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

    //color info:
    G_rgb(red[obj_index][i], grn[obj_index][i], blu[obj_index][i]) ;

    //fill em up
    G_fill_polygon(xp, yp, np) ;
  }
}



//centers object file file (800,800)
void center_object(int objnum){
  double xmin, xmax, ymin, ymax;
  double center_x, center_y;
  double dx, dy, sx, sy,nx,ny;

  //find center of bounding box string min x, min y, max x, maxy
  xmin = xmax = x[objnum][0];
  ymin = ymax = y[objnum][0];
  
  
  for(int i = 0; i < numpoints[objnum]; i++){
    if(x[objnum][i] < xmin) xmin = x[objnum][i];
    if(x[objnum][i] > xmax) xmax = x[objnum][i];
    if(y[objnum][i] < ymin) ymin = y[objnum][i];
    if(y[objnum][i] > ymax) ymax = y[objnum][i];
  }
    
  center_x = (xmin + xmax) / 2;
  center_y = (ymin + ymax) / 2;
  
  dx = (400 - center_x);
  dy = (400 - center_y);

  //matrices functions :)
  M2d_make_identity(matrix);
  M2d_make_translation(matrix, dx, dy);

  //    apply some translation 
  for (int i = 0; i < numpoints[objnum]; i++) {
    x[objnum][i] = matrix[0][0] * x[objnum][i] + matrix[0][1] * y[objnum][i] + matrix[0][2];
    y[objnum][i] = matrix[1][0] * x[objnum][i] + matrix[1][1] * y[objnum][i] + matrix[1][2];
  }
}

	   

void spin_obj(int objnum, double radians){
  double temp;

  //matrice functions :)
  M2d_make_identity(matrix);
  M2d_make_rotation_cs(matrix, cos(radians), sin(radians));

  //lets create some rotation
  for (int i = 0; i < numpoints[objnum]; i++) {
    temp = matrix[0][0] *( x[objnum][i]-400) - matrix[0][1] * (y[objnum][i]-400) + 400;
    y[objnum][i] = matrix[0][1] * (x[objnum][i]-400) + matrix[1][1] * (y[objnum][i]-400) + 400;
    x[objnum][i] = temp;
  }
}



int main(int argc, char ** argv) {
    FILE *fin;
    int key, obj_num;
    char fname[100];
    int num_objects = 0;
    double spin_angle = 4.0 *M_PI/180;

    //reading in files
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
    G_clear() ;
    G_rgb(0,0,1) ;
    draw_object(obj_num) ;

    
    //how to navigate and spin files
  int  q,k ;
  double matrix2[3][3], matrix3[3][3];
  
  while (1) {

    M2d_make_identity(matrix) ;
    
    q = G_wait_key() ;
    
    if (q == 'q') {
      exit(0) ;

    } else if (('0' <= q) && (q <= '9')) { //move through files
      printf("action: %ls\n", &q);
      k = q - '0' ;

      if(k < num_objects){obj_num = k;}

    } else if ((q == 'a') || (q == 'A')) {
      printf("action: %ls\n", &q);
      M2d_make_translation(matrix, -2, 0);


    } else if ((q == 'w') || (q == 'W')) {
      printf("action: %ls\n", &q);
      M2d_make_translation(matrix, 0, 2);

    } else if ((q == 's') || (q == 'S')) {
      printf("action: %ls\n", &q);
      M2d_make_translation(matrix, 0, -2);

    } else if ((q == 'd') || (q == 'D')){
      printf("action: %ls\n", &q);
      M2d_make_translation(matrix, 2, 0);

    } else if ((q == 'l') || (q == 'L')){
      printf("action: %ls\n", &q);
      spin_obj(obj_num, spin_angle);

      
    }
    else if ((q == 'r') || (q == 'R')){
      printf("action: %ls\n", &q);
      M2d_make_translation(matrix,  -x[obj_num][numpoints[obj_num]], -y[obj_num][numpoints[obj_num]]);
      M2d_make_rotation_cs(matrix2, sin(spin_angle), cos(spin_angle));
      M2d_mat_mult(matrix2, matrix2, matrix);
      M2d_make_translation(matrix3,  x[obj_num][numpoints[obj_num]], y[obj_num][numpoints[obj_num]]);
      M2d_mat_mult(matrix, matrix3, matrix2);


    }else if ((q == 'R') || (q == 'R')){
      printf("action: %ls\n", &q);
      M2d_make_translation(matrix, 2, 0);

    }else {
      printf("no action\n") ;
    }

    M2d_mat_mult_points (x[obj_num], y[obj_num],  matrix,
			 x[obj_num], y[obj_num], numpoints[obj_num]+1) ;

    G_rgb(0,0,0) ; 
    G_clear() ;
    G_rgb(0,0,1) ;
    draw_object(obj_num) ;
    
  }
}
  
