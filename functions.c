#include "FPToolkit.c"
  
void my_polygon(double x[3], double y[3], int n){
  // G_rgb(0,0,0); //black
}
  
  
int click_and_save(double x2[3], double y2[3]){
    double p[2];
    G_rgb(1,0,0); G_fill_rectangle(0,0,20,10);
    int n = 0;

    while(n <= 10){
      G_wait_click(p);
      G_circle(p[0], p[1], 2);
      G_line(x2[n], y2[n], x2[n+1], y2[n+1]);
      x2[n] = p[0]; y2[n] = p[1];
      n++; }

    return n;
 }


int main(){
  double a[100], b[100], c[100], d[100];
  double ncd, nab;
  int swidth = 800;
  int sheight = 800;

  G_init_graphics(swidth,sheight);

  G_rgb(0.3,0.3,0.3);
  G_clear();
  
  //nab = my_polygon(a, b, 4);
  ncd = click_and_save(c,d);
  int G_wait_click(double v[2])
  nab = click_and_save(a,b);

  my_polygon(a, b, nab);
  my_polygon(c, d, ncd);

}
  