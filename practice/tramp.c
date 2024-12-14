
#include "FPToolkit.c"
#include <stdbool.h>
#define MAXPTS 59000
#define MAXPOLYS 57500


int numpoints ;
int numpolys ;
double x[MAXPTS] ;
double y[MAXPTS] ;
int psize[MAXPOLYS] ;
int con[MAXPOLYS][20] ;
double red[MAXPOLYS],grn[MAXPOLYS],blu[MAXPOLYS] ;



int read_object(FILE *f){

    // reads number of points from a file
    fscanf(f,"%d",&numpoints) ;

    //if the point info is more than maxpoints
    if (numpoints >= MAXPTS) {
      // need an extra for object centering
      printf("MAXPTS = %d :  exceeded.\n",MAXPTS) ;
      //kill
      exit(1) ;
    }

    for (int i = 0 ; i < numpoints ; i++) {
      //scans to see the x and y coordinates of the file
      fscanf(f,"%lf %lf",&x[i],&y[i]) ;
    }

    // connectivity info
    //reads number of polygons
    fscanf(f,"%d",&numpolys) ;
    //check to see if n of polygons is higher than MaxPoly
    if (numpolys > MAXPOLYS) {
      printf("MAXPOLYS = %d :  exceeded.\n",MAXPOLYS) ;
      //if so kill
      exit(1) ;
    }
    
    
    for (int i = 0 ; i < numpolys ; i++) {
      //scan for how many vertices the polygon has
      fscanf(f,"%d",&psize[i]) ;
      for (int j = 0 ; j < psize[i] ; j++) {
	//checking for the storage of indices of vertices
	//that make up each polygon in a 2d array
        fscanf(f,"%d",&con[i][j]) ;
      } 
    } 

    
    // color info :
    for (int i = 0 ; i < numpolys ; i++) {
      fscanf(f,"%lf %lf %lf",&red[i],&grn[i],&blu[i]) ;
    }    
}

void center_object(int objnum){
  //find center of bounding box sing min x, min y, max x, maxy
  //transiton center poitn to orgin (0,0)
  //then to center screen (400,400)
  //apply some transiont maax to all x[], and y[]
  //in the shape to move whole thing
}


//void translate(int onum, double dx, double dy){
//int numpoints[2];
//int x[2][2];
//int y[2][2];
//for(int i = 0; i < numpoints[onum]; i++){
//  x[onum][i] += dx;
//  y[onum][i] += dy;
//}
//}

//void scale(int onum, double sx, double sy){
// int numpoints[2];
//int x[2][2];
//int y[2][2];
//for(int i = 0; i < numpoints[onum]; i++){
//  x[onum][i] *= sx;
//  y[onum][i] *= sy;
//}
//}
	    


//void spin_obj(int onum, double degrees){
//int i;
//int numpoints[2];
//int x[2][2];
//int y[2][2];
//double t, temp, c, s;
//c = cos(t);
//s = sin(t);

//for(i = 0; i < numpoints[onum]; i++){
//  temp = x[onum][i]*c - y[onum][i]*s;
//  y[onum][i] = x[onum][i]*s + y[onum][i]*c;
//  x[onum][i] = temp;
//}
//}


int print_object (FILE *fout){

  fprintf(fout, "%d\n",numpoints) ;

  for (int i = 0 ; i < numpoints ; i++) {
    //12.6 ensures things are printed out with a
    //width of 12 charactesr includeing 6 decimal
    fprintf(fout, "%12.6lf %12.6lf\n",x[i],y[i]) ;
  }
  

  for (int i = 0 ; i < numpolys ; i++) {
    //integer shoud occupy at least three spaces
    fprintf(fout, "%3d    ",psize[i]) ;

    for (int j = 0 ; j < psize[i] ; j++) {
      //integer shoudl occupy at least two spaces
      fprintf(fout, "%2d ", con[i][j]) ; }
    
    fprintf(fout, "\n") ;
  }


  for (int i = 0 ; i < numpolys ; i++) {
    fprintf(fout,"%lf %lf %lf\n",red[i],grn[i],blu[i]) ;
  }      
}


int draw_object ()
{
  int h;
  double xp[100],yp[100] ;
  int np ;

  for (int i = 0 ; i < numpolys ; i++) {

    np = psize[i] ; //get number of vertices in this polygon
    
    //get vertix coordinates for current polygon
    for (int j = 0 ; j < np ; j++) {
      h = con[i][j] ; // get index of vertex
      xp[j] = x[h] ; //store the x coordinate
      yp[j] = y[h] ; //store y coordiante
    }

    G_rgb(red[i], grn[i], blu[i]) ;
    G_fill_polygon(xp,yp,np) ;
  }
}


int main(){
  
  FILE *fin ;
  int key, objnum = -1;
  char fname[100] ;
  double click_point[2];
  

  printf("enter names of xy file ") ;
  scanf("%s",fname) ;

  //open file to read 'r'
  fin = fopen(fname,"r") ;
  
  if (fin == NULL) {
      printf("can't read file, %s\n",fname) ;
      exit(1) ;
  }

  read_object(fin) ;


  //  print_object(stdout, 0) ;
  G_init_graphics(800, 800);
  

  while(true){
    G_rgb(0,0,0) ;
    G_clear() ;
    draw_object() ;
    key = G_wait_click(click_point);

    
    if (key >= 'a' && key <= 'o') {
      int temp = key - 'a';
	
      if(temp != objnum){
	objnum = temp;
	G_rgb(0,0,0) ;
	G_clear() ;
	draw_object() ;
        key = G_wait_click(click_point);     }
    }
  
    if(key == 'g'){
      break;
    }
  }

  // int numobjects;
  
  // for(int k = 0; k < numobjects; k++){
  //center_object(k);
    //some function to initilize object size (kinda optional)
  //}

  //fclose(fin);
  //return 0;

}
