//////////////////////////////////////////////////////////////////////////////
////////////////////////Phong reflection model//////////////////////////////////
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
double color[3][3] = {{1,1,1},
		 {0,1,0},
		 {0,0,1}};

double matrix[4][4] = {0};
THING listall[MAXPTS];
double lco[3] = {800.0, 800.0, 0.0};
double ambient = 0.4, diffuseMax = 0.5, specularPow = 35;

int film (double xp[], double yp[], double zp[], int pts, double flmX[], double flmY[]) {
  double xTemp[pts], yTemp[pts] ;

  for (int i = 0; i < pts; i++) {
    flmX[i] = ((xp[i]/zp[i]) * 800/2) + 800/2 ;
    flmY[i] = ((yp[i]/zp[i]) * 800/2) + 800/2 ;
  }
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



double dot_product(double a[], double b[])
/* Own version of a dot product */
{
  double result = 0;

  for(int i = 0; i < 3; i++){
    result = result + a[i] * b[i];
  }

  return result;

}

void normalize_vector(double vec[3]) {
  /* normalizes vecotrs */
  double magnitude = sqrt(dot_product(vec, vec));
  if (magnitude != 0) {
    vec[0] /= magnitude;
    vec[1] /= magnitude;
    vec[2] /= magnitude;
    }
}

void normal_v(double xx[], double yy[], double zz[], double normal[3])
/* function to find normal vector */
{
  double E[3], v1[3], v2[3];

  v1[0] = xx[1] - xx[0];
  v1[1] = yy[1] - yy[0];
  v1[2] = zz[1] - zz[0];

  v2[0] = xx[2] - xx[0];
  v2[1] = yy[2] - yy[0];
  v2[2] = zz[2] - zz[0];

  //calculate normal vector
  M3d_x_product(normal, v1, v2);
  normalize_vector(normal);
}


void view_v(double xx[], double yy[], double zz[], double view[3])
/* function to find viewing vector */
{
  double centroid[3];
  //calculate the center
  centroid[0] = (xx[0] + xx[1] + xx[2])/3;
  centroid[1] = (yy[0] + yy[1] + yy[2])/3;
  centroid[2] = (zz[0] + zz[1] + zz[2])/3;
  //calculate the center
  view[0] = 0 - centroid[0];
  view[1] = 0 - centroid[1];
  view[2] = 0 - centroid[2];
  normalize_vector(view);
}


double light_v(double xx[], double yy[], double zz[], double ka[])
/* function to find light vector */
{
  double centroid[3];
  //calculate the center
  centroid[0] = (xx[0] + xx[1] + xx[2])/3;
  centroid[1] = (yy[0] + yy[1] + yy[2])/3;
  centroid[2] = (zz[0] + zz[1] + zz[2])/3;

  ka[0] = ka[0] + centroid[0];
  ka[1] = ka[1] + centroid[1];
  ka[2] = ka[2] + centroid[2];

  normalize_vector(ka);
}



double light_em_up(double xx[], double yy[], double zz[], int pnum,
		 double ka, double kd, double alpha)
/*
   pnum: number of points in a polygon
   ka: ambient reflection constant
   kd: diffuse reflection constant
   alpha: shineness constant
*/

{
  double normal[3], view[3], light[3], reflection[3];
  normal_v(xx, yy, zz, normal);
  view_v(xx, yy, zz, view);
  light_v(xx, yy, zz, light);

  double n_dot_l = dot_product(normal, light);
  reflection[0] = 2.0 * n_dot_l * normal[0] - light[0];
  reflection[1] = 2.0 * n_dot_l * normal[1] - light[1];
  reflection[2] = 2.0 * n_dot_l * normal[2] - light[2];

  normalize_vector(reflection);


  double ambient = ka;
  double diffuse_term = kd  * (dot_product(normal, light));

  double specular_term = (1 - (ambient + kd)) *
    (pow(dot_product(view, reflection), alpha));

  double total_intensity = ambient + diffuse_term + specular_term;

  // Debugging print statements
  printf("Light Vector: (%f, %f, %f)\n", light[0], light[1], light[2]);
  printf("View Vector: (%f, %f, %f)\n", view[0], view[1], view[2]);
  printf("Normal Vector: (%f, %f, %f)\n", normal[0], normal[1], normal[2]);
  printf("Ambient Term: %f\n", ambient);
  printf("Diffuse Term: %f\n", diffuse_term);
  printf("Refelction Vector: (%f, %f, %f)\n", reflection[0], reflection[1], reflection[2]);
  printf("Specular Term: %f\n", specular_term);
  printf("Total Intensity: %f\n", total_intensity);

  if(ambient > total_intensity){return ambient;}
  else{return total_intensity;}


}


int backface_removal(double xx[], double yy[], double zz[], int n)
{
  double v1[3]; //edge vector
  double v2[3]; //edge vector
  double normal[3]; //normal vector
  double sol; //var for storing dot product
  double E[3]; //viewing vector
  int signal; //whether the vector should or should not be drawn ;)


  //edge vectors
  v1[0] = xx[1] - xx[0];
  v1[1] = yy[1] - yy[0];
  v1[2] = zz[1] - zz[0];

  v2[0] = xx[2] - xx[0];
  v2[1] = yy[2] - yy[0];
  v2[2] = zz[2] - zz[0];

  //calculate normal vector
  M3d_x_product(normal, v1, v2);

  //viewing vector
  E[0] = 0 - xx[0];
  E[1] = 0 - yy[0];
  E[2] = 0 - zz[0];

  sol = dot_product(normal, E);
  printf("%d \n", rvs);
  //allows user to invert from front-face to back-face
  sol = (sol * rvs);

  //if more than 0 draw lines...
  if(sol > 0){signal = 1;}

  //dont draw
  else{signal = 0;}

  //return so you know
  return signal;
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

  x[objnum][numpoints[objnum]] = 0;
  y[objnum][numpoints[objnum]] = 0;
  z[objnum][numpoints[objnum]] = 0;
}



// draws an individual polygon from an object
int draw_single_poly (int objnum, int polynum) {
  int np, h, i ;
  double xp[100], yp[100], zp[100] ;
  double xFilm[100], yFilm[100] ;
  double color1, color2, color3, flux, intensity, signal = 0 ;

  np = psize[objnum][polynum] ;

  for (i = 0; i < np; i++) {
    h = con[objnum][polynum][i] ;
    xp[i] = x[objnum][h] ;
    yp[i] = y[objnum][h] ;
    zp[i] = z[objnum][h] ;
  }

 
  intensity = light_em_up(xp, yp, zp, np, ambient, diffuseMax, specularPow) ;
  
  color1 = color[objnum][0] ;
  color2 = color[objnum][1] ;
  color3 = color[objnum][2] ;
  flux = ambient + diffuseMax ;

  if (intensity <= flux) {
    flux = intensity / flux ;
    color1 = color[objnum][0] * flux ;
    color2 = color[objnum][1] * flux ;
    color3 = color[objnum][2] * flux ;
    
  } else {
    flux = (intensity - flux) / (1 - flux) ;
    color1 = color[objnum][0] + flux*(1 - color[objnum][0]) ;
    color2 = color[objnum][1] + flux*(1 - color[objnum][1]) ;
    color3 = color[objnum][2] + flux*(1 - color[objnum][2]) ;
  }

  film(xp, yp, zp, np, xFilm, yFilm) ;
  signal = backface_removal(xp, yp, zp, np);

    //if 1 draw it
  if(signal == 1){
    G_rgb(color1,color2,color3) ;
    G_fill_polygon(xFilm, yFilm, np) ;
    G_rgb(0, 0, 0) ;
    G_polygon(xFilm, yFilm, np) ;
  }
}



int init_array(int nah)
{
  int np = 0;
  int indx = 0;
  int v = 0;
  double distance = 0;

  for(int onum = 0; onum < nah; onum++){
    for(int pnum = 0; pnum < numpolys[onum]; pnum++){
      np = psize[onum][pnum];

      listall[indx].objnum = onum;
      listall[indx].polynum = pnum;

      for (int j = 0; j < np; j++) {
	v = con[onum][pnum][j];
	distance += sqrt(pow(x[onum][v], 2) + pow(y[onum][v], 2)+
			 pow(z[onum][v], 2));
      }

      listall[indx].dist = distance / np;
      indx++;
      distance = 0;

    }
  }
  return indx;
}




void print_array(int indx)
{
  int i ;
  for (i = 0 ; i < indx ; i++) {
    printf("objnum:%d polynum:%d dist:%lf\n",listall[i].objnum,
	   listall[i].polynum, listall[i].dist) ;
  }
  printf("\n");
}



static int compare(const void *p, const void *q)
{
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return 1 ;
  else if (((*a).dist) > ((*b).dist)) return -1 ;
  else return 0 ;

}


void draw_all_object(int nah)
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
