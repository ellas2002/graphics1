#include "FPToolkit.c"

void intersection(double a[], double b[], double c[]){
  double mxy;
  double maxy;
  double d[2];

  if(a[1] < b[1]){ mxy = a[1]; maxy = b[1];}
  
  else{mxy = b[1]; maxy = a[1];}

  if(mxy < c[1] && c[1] < maxy){
  d[1] = c[1];

  d[0] = a[0] + (c[1] - a[1]) * (b[0] - a[0]) /(b[1] - a[1]);

  
  G_fill_circle(d[0], d[1], 3);}

  else{printf("does not intersect");}

}


int main(){

  G_init_graphics(800,800);
  G_rgb(2,1,2);
  G_clear();

  G_rgb(1,0,0);
  
  double x[2], y[2];
  
  G_wait_click(x);
  G_fill_circle(x[0],x[1],3);

  G_wait_click(y);
  G_fill_circle(y[0],y[1],3);


  G_line(x[0],x[1],y[0],y[1]);

  double p[2];
  
  G_wait_click(p);
  G_fill_circle(p[0], p[1], 3);

  G_line(0, p[1], 800, p[1]);

  double i[2];
  double mxy, maxy;

  intersection(x, y, p);
       
  int key;
  key = G_wait_key();

}
