//////////////////////////////////////////////////////////////////////////////
////////////////////////back-face elimination/////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include <stdbool.h>
#include <math.h>
#define MAXPTS 59000
#define MAXPOLYS 57500
#define MAX_OBJECTS 10


int numpoints[MAX_OBJECTS] ;
int numpolys[MAX_OBJECTS] ;
double x[MAX_OBJECTS][MAXPTS] ;
double y[MAX_OBJECTS][MAXPTS] ;
double z[MAX_OBJECTS][MAXPTS] ;
int psize[MAX_OBJECTS][MAXPOLYS] ;
int con[MAX_OBJECTS][MAXPOLYS][20] ;
double red[MAX_OBJECTS][MAXPOLYS],grn[MAX_OBJECTS][MAXPOLYS],blu[MAX_OBJECTS][MAXPOLYS] ;
double matrix[4][4] = {0};


//reads the file in order to extract info and output the file
int read_object(int obj_index, FILE *f)
{

  //scan in number of points in object
  fscanf(f,"%d",&numpoints[obj_index]) ;

  //if the number of points is more than the alloted amount: kill
  if (numpoints[obj_index] >= MAXPTS) {
      printf("MAXPTS = %d :  exceeded.\n",MAXPTS) ;
      exit(1) ;
  }

  //scan in the x and y points
  for (int i = 0 ; i < numpoints[obj_index] ; i++) {
    fscanf(f,"%lf %lf %lf",&x[obj_index][i],&y[obj_index][i], &z[obj_index][i]) ;
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
      fscanf(f,"%lf",&red[obj_index][i]) ;
   }    
}




double dot_product(double a[], double b[])
{
  double result = 0;
  
  for(int i = 0; i < 3; i++){
    result = result + a[i] * b[i];
  }
  
  return result;

}



int backface_removal(double xx[], double yy[], double zz[], int n)
{
  double v1[3]; //edge vector
  double v2[3]; //edge vector
  double normal[3]; //normal vector
  double sol; //var for storing dot product
  double E[3]; 
  int rvs = -1; //reversal factor
  int signal; //whether the vector should or should not be drawn ;)

  E[0] = 0 - xx[0];
  E[1] = 0 - yy[0];
  E[2] = 0 - zz[0];

  //find edge vectors 
  v1[0] = xx[1] - xx[0];
  v1[1] = yy[1] - yy[0];
  v1[2] = zz[1] - zz[0];
  
  v2[0] = xx[2] - xx[0];
  v2[1] = yy[2] - yy[0];
  v2[2] = zz[2] - zz[0];

  //calculate normal vector
  M3d_x_product(normal, v1, v2);
  
  sol = dot_product(normal, E);
  sol = (sol * rvs);

  //if more than 0 draw lines...
  if(sol > 0){signal = 1;}

  //dont draw
  else{signal = 0;}

  //return so you know
  return signal;  
}



//draws the object in the scanned file
int draw_object(int obj_index)
{
  double xpp[100], ypp[100], zpp[100], xp[100], yp[100], zp[100], H, h ;
  int np ;
  int signal;
  h = 40.0 * M_PI/180;
  
  H  = tan(h);
  
  for (int i = 0 ; i < numpolys[obj_index] ; i++) {
    np = psize[obj_index][i] ; //get number of vertices in this polygon
    
    //get vertix coordinates for current polygon
    for (int j = 0 ; j < np ; j++) {
      int n = con[obj_index][i][j] ; // get index of vertex
      xpp[j] = (x[obj_index][n]/ z[obj_index][n])  * (400/H) + 400  ; //store the x coordinate
      ypp[j] = (y[obj_index][n]/ z[obj_index][n]) * (400/H) + 400 ; //store y coordiante
      
      xp[j] = x[obj_index][n];
      yp[j] = y[obj_index][n];
      zp[j] = z[obj_index][n];
    }
    //find the vectors thats either 0 or 1
    signal = backface_removal(xp, yp, zp, np);

    //if 1 draw it
    if(signal == 1){
      G_polygon(xpp, ypp, np) ;
    }
  }
}


//centers object file file (800,800)
void center_object(int objnum)
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double center_x, center_y, center_z;
  double dx, dy,dz, sx, sy,nx,ny;

  //find center of bounding box string min x, min y, max x, maxy
  xmin = xmax = x[objnum][0];
  ymin = ymax = y[objnum][0];
  zmin = zmax = z[objnum][0];
  
  
  for(int i = 0; i < numpoints[objnum]; i++){
    if(x[objnum][i] < xmin) xmin = x[objnum][i];
    if(x[objnum][i] > xmax) xmax = x[objnum][i];
    if(y[objnum][i] < ymin) ymin = y[objnum][i];
    if(y[objnum][i] > ymax) ymax = y[objnum][i];
    if(z[objnum][i] < zmin) zmin = z[objnum][i];
    if(z[objnum][i] > zmax) zmax = z[objnum][i];
  }
    
  center_x = (xmax - xmin) / 2;
  center_y = (ymax - ymin) / 2;
  center_z = (zmax - zmin) / 2;
  
  dx = center_x + xmin;
  dy = center_y + ymin;
  dz = center_z + zmin;

  //matrices functions :)

  M3d_make_identity(matrix);
  M3d_make_translation(matrix, dx, dy, dz);
      
  x[objnum][numpoints[objnum]] = 0;
  y[objnum][numpoints[objnum]] = 0;
  z[objnum][numpoints[objnum]] = 0;
}
	



int main(int numFiles, char ** file)
{
  FILE *fin;
  int key;
  char fname[100];
  int num_objects = 0;
  double spin_angle = 15.0 * M_PI/180;
  double matrix2[4][4], matrix3[4][4];
  int  sign = 1 , action = 't' ;  
  int  onum = 0 ;
  int rvs = -1;
  int  q, k ;

  //reading in files
  for(int i = 0; i < numFiles; i++){
    if(i + 1 < numFiles){
      printf("argv[%d] is %s\n", i + 1, file[i + 1]);
    }
  }

  num_objects = numFiles -  1;
    
  //open files (read)
  for(onum = 1; onum < num_objects + 1; onum++){
    fin = fopen(file[onum], "r") ;

    //if no files
    if (fin == NULL) {
      printf("can't read file, %s\n", file[onum]) ;
      exit(1) ;
    }
    
    read_object(onum - 1, fin);
    center_object(onum);
    fclose(fin);
  }

  //initialize
  G_init_graphics(800, 800);
  G_clear();
  G_rgb(1,0,0);
  draw_object(0);
  onum = 0;

  while (1) {
    M3d_make_identity(matrix);

    q = G_wait_key() ;      
    
    if (q == 'q') { //exit
      exit(0) ;

    } else if (q == 'c') { //flips direction
      printf("action: %ls\n", &q);
      sign = -sign ;

    } else if (q == 'i') { //inverts backface culling
      printf("action: %ls\n", &q);
      rvs = -1 * rvs ;

    } else if (q == 't') { //for translating
      printf("action: %ls\n", &q);
      action = q ;

    } else if (q == 'r') { //for rotating
      printf("action: %ls\n", &q);
      action = q ;

    } else if (('0' <= q) && (q <= '9')) { //moving between files
      k = q - '0' ;  
      if (k < num_objects) { onum = k ; }

    } else if ((q == 'x') && (action == 't')) { // moves left or right
      printf("action: %ls\n", &q);
      M3d_make_translation(matrix, sign * 0.2, 0, 0);

    } else if ((q == 'y') && (action == 't')) { //moves up or down
      printf("action: %ls\n", &q);
      M3d_make_translation(matrix, 0, sign * 0.2, 0);

    } else if ((q == 'z') && (action == 't')) { //moves forward or backwards
      printf("action: %ls\n", &q);
      M3d_make_translation(matrix, 0, 0, sign * 0.2);

    } else if ((q == 'x') && (action == 'r')) { //rotates on x axis
      printf("action: %ls\n", &q);
      M3d_make_translation(matrix,  -x[onum][numpoints[onum]], -y[onum][numpoints[onum]], -z[onum][numpoints[onum]]);
      M3d_make_x_rotation_cs(matrix2, cos(spin_angle * sign), sin(spin_angle * sign));
      M3d_mat_mult(matrix2, matrix2, matrix);
      M3d_make_translation(matrix3,  x[onum][numpoints[onum]], y[onum][numpoints[onum]], z[onum][numpoints[onum]]);
      M3d_mat_mult(matrix, matrix3, matrix2);
      

      
    } else if ((q == 'y') && (action == 'r')) { // rotates on y axis
      printf("action: %ls\n", &q);
      M3d_make_translation(matrix,  -x[onum][numpoints[onum]], -y[onum][numpoints[onum]], -z[onum][numpoints[onum]]);
      M3d_make_y_rotation_cs(matrix2, cos(spin_angle*sign), sin(spin_angle*sign));
      M3d_mat_mult(matrix2, matrix2, matrix);
      M3d_make_translation(matrix3,  x[onum][numpoints[onum]], y[onum][numpoints[onum]], z[onum][numpoints[onum]]);
      M3d_mat_mult(matrix, matrix3, matrix2);
      
    } else if ((q == 'z') && (action == 'r')) { // rotates on z axis
      printf("action: %ls\n", &q);
      M3d_make_translation(matrix,  -x[onum][numpoints[onum]], -y[onum][numpoints[onum]], -z[onum][numpoints[onum]]);
      M3d_make_z_rotation_cs(matrix2, cos(spin_angle*sign), sin(spin_angle*sign));
      M3d_mat_mult(matrix2, matrix2, matrix);
      M3d_make_translation(matrix3,  x[onum][numpoints[onum]], y[onum][numpoints[onum]], z[onum][numpoints[onum]]);
      M3d_mat_mult(matrix, matrix3, matrix2);

      
    } else {printf("no action\n") ;}

    
    M3d_mat_mult_points (x[onum],y[onum],z[onum],  matrix,
			 x[onum],y[onum],z[onum],numpoints[onum]+1) ;

    G_rgb(0,0,0) ; 
    G_clear() ;
    G_rgb(1,0,0) ;
    
    
    draw_object(onum);
  }
}
