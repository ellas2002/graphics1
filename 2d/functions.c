#include "FPToolkit.c"
#include <stdio.h>

//swap moves one number to another spot
void swap(double *a, double *b){
  int temp = *a;
  *a = *b;
  *b = temp;
}

//sort chages the posiion of a number depending on how large or small it is, to put into order
void sort(double v[], int n) {
    for (int i = 0; i < n - 1; i++) {
        int min = i;
        
        for (int j = i + 1; j < n; j++) {
            if (v[j] < v[min]) {
                min = j;
            }
        }
        
        if (min != i) {
            swap(&v[min], &v[i]);
        }
    }
}

//printArray allows me to see if the sorting algo works
void printArray(double arr[], int size)
{
    int i;
    
    for (i=0; i < size; i++){printf("%f ", arr[i]);}
    
    printf("\n");
}


//intersection finds the point of intersection between two lines 
double intersection(double a[], double b[], double c[], double n[]){
  double mxy;
  double maxy;
  double d[2];

  if(a[1] < b[1]){ mxy = a[1]; maxy = b[1];}
  
  else{mxy = b[1]; maxy = a[1];}

  if(mxy < c[1] && c[1] < maxy){
  d[1] = c[1];

  d[0] = a[0] + (c[1] - a[1]) * (b[0] - a[0]) /(b[1] - a[1]);

  
  G_fill_circle(d[0], d[1], 3);}

  //else{printf("does not intersect");}

  return 1;

}

// saves the clicks the user inputs  
int click_and_save(double x[], double y[]){
  G_rgb(0,1,0);
  G_fill_rectangle(0,0,20,10);
  
  double v[2];
  int n;
  n = 0;

  while(0 == 0){
      G_wait_click(v);
      
      if( v[0] < 20 && v[1] < 10) {break;}

      x[n] = v[0]; y[n] = v[1]; n++;

      G_fill_circle(v[0],v[1],2);
    }
  
  return n;
 }


void my_polygon(double x[], double y[], int n){
  int i, j;
  i = 0;
  
  while( i < n){
    j = i + 1;
    
    if(j == n){j=0;}
    
    G_line(x[i], y[i], x[j], y[j]);
    
    i++;
  }
  G_wait_key();

}

void my_fill_polygon(double x[], double y[], int n){
    // Iterate through each scanline
    for (double ycheck = 0.1; ycheck < 800; ycheck += 0.5) {
        double intersections[100];
        int counter = 0;

        // Find intersections with all polygon edges
        for(int i = 0; i < n; i++) {
	  //navigates through the vertices of the polygon in a circular manner
            int j = (i + 1) % n;
	    //if y[i] is in range of ycheck
	    if ((y[i] <= ycheck && y[j] > ycheck) || (y[j] <= ycheck && y[i] > ycheck)) {
                // Calculate x-coordinate of intersection using slope formula
                double x_intersect = x[i] + (ycheck - y[i]) * (x[j] - x[i]) / (y[j] - y[i]);
                intersections[counter++] = x_intersect;
            }
        }

        // Sort the intersections
	sort(intersections, counter);

        // Draw horizontal lines between pairs of intersections
	for (int i = 0; i < counter; i += 2) {
	  G_line(intersections[i], ycheck, intersections[i+1], ycheck);
	  //fill em up
	  G_wait_key();	
	}
    }	  	
}



int main(){
  double a[100], b[100];
  int nab;
  
  G_init_graphics(800, 800); //initializes
  G_rgb(0.3,0.3,0.3);
  G_clear();
  
  G_rgb(0,1,0);
  
  nab = click_and_save(a, b);
  my_polygon(a, b , nab);

  my_fill_polygon(a, b, nab);
  
  G_wait_key();

}
    
  

