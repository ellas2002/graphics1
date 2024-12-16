//////////////////////////////////////////////////////////////////////////////
////////////////////////painters algorithm///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "../FPToolkit.c"
#include "M3d_matrix_tools.c"
#define MAXPTS 59000
#define MAXPOLYS 57500
#define MAX_OBJECTS 10

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}

  
THING;
int numpoints[MAX_OBJECTS] ; //points
int numpolys[MAX_OBJECTS] ; //polygon
int rvs = 1; //reversal factor
double x[MAX_OBJECTS][MAXPTS] ;
double y[MAX_OBJECTS][MAXPTS] ;
double z[MAX_OBJECTS][MAXPTS] ;
int psize[MAX_OBJECTS][MAXPOLYS] ; //polygon vertices
int con[MAX_OBJECTS][MAXPOLYS][20] ; //indices of vertices
double red[MAX_OBJECTS][MAXPOLYS],grn[MAX_OBJECTS][MAXPOLYS],blu[MAX_OBJECTS][MAXPOLYS] ;
double color[3][3] = {{1,0,0},
		 {0,1,0},
		 {0,0,1}};

double matrix[4][4] = {0};
THING listall[MAXPTS];


int film (double xp[], double yp[], double zp[], int pts, double flmX[], double flmY[]) {
  double xTemp[pts], yTemp[pts] ;

  for (int i = 0; i < pts; i++) {
    flmX[i] = ((xp[i]/zp[i]) * 800/2) + 800/2 ;
    flmY[i] = ((yp[i]/zp[i]) * 800/2) + 800/2 ;
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

  M3d_make_translation(matrix,  x[objnum][numpoints[objnum]],
		       y[objnum][numpoints[objnum]],
		       z[objnum][numpoints[objnum]]);

  x[objnum][numpoints[objnum]] = 0;
  y[objnum][numpoints[objnum]] = 0;
  z[objnum][numpoints[objnum]] = 0;
}



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
}



// draws an individual polygon from an object
int draw_single_poly (int objnum, int polynum) {
  /* draws a single polygon function, that takes
     into acount the film and color of said polygon */
  int np, h, i ;
  double xp[100], yp[100], zp[100] ;
  double xFilm[100], yFilm[100] ;
  
  np = psize[objnum][polynum] ;

  for (i = 0; i < np; i++) {
    h = con[objnum][polynum][i] ;
    xp[i] = x[objnum][h] ;
    yp[i] = y[objnum][h] ;
    zp[i] = z[objnum][h] ;
  }

  film(xp, yp, zp, np, xFilm, yFilm) ;
  
  G_rgb(color[objnum][0],color[objnum][1],color[objnum][2]) ;
  G_fill_polygon(xFilm, yFilm, np) ;
  G_rgb(0, 0, 0) ;
  G_polygon(xFilm, yFilm, np) ;
}


int init_array(int nah)
/* sets up THING with actual scanned in paramaters,
   measures the depth of z coordniate  */
{
  int np = 0;
  int indx = 0;
  int v = 0;
  double ave_z = 0;
  

  for(int onum = 0; onum < nah; onum++){
    for(int pnum = 0; pnum < numpolys[onum]; pnum++){
      np = psize[onum][pnum];

      listall[indx].objnum = onum;
      listall[indx].polynum = pnum;
      
      for (int j = 0; j < np; j++) {
	v = con[onum][pnum][j];
	ave_z += z[onum][v];
      }
      
      listall[indx].dist = ave_z / np;
      indx++;
      ave_z = 0;
      
    }
  }
  return indx;
}

  


void print_array(int indx)
/* prints the init_array */
{
  int i ;
  for (i = 0 ; i < indx ; i++) {
    printf("objnum:%d polynum:%d dist:%lf\n",listall[i].objnum,
	   listall[i].polynum, listall[i].dist) ;
  }
  printf("\n");
}


  
static int compare(const void *p, const void *q)
/* compare function for depth of z */
{    
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return 1 ;
  else if (((*a).dist) > ((*b).dist)) return -1 ;
  else return 0 ;

}


void draw_all_object(int nah)
/* sorts all the information in listall and
   compares them based on category.
 prints out the the polygons with this info*/
{
    int total_polygons = init_array(nah);
    
    qsort(listall, total_polygons, sizeof(THING), compare);
    
  
    for (int i = 0; i < total_polygons; i++) {
      draw_single_poly(listall[i].objnum, listall[i].polynum);
    }
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
  onum = 0;

  while (1) {
    M3d_make_identity(matrix);

    q = G_wait_key() ;      
    
    if (q == 'q') { //exit
      exit(0) ;

    } else if (q == 'c') { //flips direction
      sign = -sign ;

    } else if (q == 'i') { //inverts backface culling
      rvs = -rvs ;

    } else if (q == 't') { //for translating
      action = q ;

    } else if (q == 'r') { //for rotating
      action = q ;

    } else if (('0' <= q) && (q <= '9')) { //moving between files
      k = q - '0' ;  
      if (k < num_objects) { onum = k ; }

    } else if ((q == 'x') && (action == 't')) { // moves left or right
      M3d_make_translation(matrix, sign * 0.2, 0, 0);

    } else if ((q == 'y') && (action == 't')) { //moves up or down
      M3d_make_translation(matrix, 0, sign * 0.2, 0);

    } else if ((q == 'z') && (action == 't')) { //moves forward or backwards
      M3d_make_translation(matrix, 0, 0, sign * 0.2);

    } else if ((q == 'x') && (action == 'r')) { //rotates on x axis
      M3d_make_translation(matrix,  -x[onum][numpoints[onum]], -y[onum][numpoints[onum]], -z[onum][numpoints[onum]]);
      M3d_make_x_rotation_cs(matrix2, cos(spin_angle * sign), sin(spin_angle * sign));
      M3d_mat_mult(matrix2, matrix2, matrix);
      M3d_make_translation(matrix3,  x[onum][numpoints[onum]], y[onum][numpoints[onum]], z[onum][numpoints[onum]]);
      M3d_mat_mult(matrix, matrix3, matrix2);
      

      
    } else if ((q == 'y') && (action == 'r')) { // rotates on y axis
      M3d_make_translation(matrix,  -x[onum][numpoints[onum]], -y[onum][numpoints[onum]], -z[onum][numpoints[onum]]);
      M3d_make_y_rotation_cs(matrix2, cos(spin_angle*sign), sin(spin_angle*sign));
      M3d_mat_mult(matrix2, matrix2, matrix);
      M3d_make_translation(matrix3,  x[onum][numpoints[onum]], y[onum][numpoints[onum]], z[onum][numpoints[onum]]);
      M3d_mat_mult(matrix, matrix3, matrix2);
      
    } else if ((q == 'z') && (action == 'r')) { // rotates on z axis
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
    
    
    draw_all_object(num_objects);
  
  }
}



