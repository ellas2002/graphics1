#include "FPToolkit.c"
#include "M2d_matrix_tools.c"


int main()
{
  double a[2], b[2], dx, dy, scale, cs, sn, len;
  
  //blank matrix
  double m[3][3], m1[3][3], m2[3][3], m3[3][3], m4[3][3];
  
  // rocket
  double rx[8] = {0, 16,  7,  7,  0, -7, -7, -16 } ;
  double ry[8] = {0,  0, 15, 35, 50, 35, 15,   0 } ;

  G_init_graphics(700,700) ;  

  //initialization
  G_rgb(0,0,0) ;
  G_clear() ;
  G_rgb(0,1,1) ;
  G_fill_polygon(rx,ry,8) ;

  //input a line that will be the initial coordinates of rocket
  G_rgb(0,1,0) ;

  G_wait_click(a) ;
  G_fill_circle(a[0],a[1],2) ;

  G_wait_click(b) ;
  G_fill_circle(b[0],b[1],2) ;
  G_line(a[0], a[1], b[0], b[1]);

  dx = b[0] - a[0];
  dy = b[1] - a[1];
  len = sqrt(dx * dx + dy * dy);
  scale = len/50;
  cs = dx/len;
  sn = dy/len;

  
  M2d_make_scaling(m1,scale,scale);
  M2d_make_rotation(m2, -M_PI/2);
  M2d_make_rotation_cs(m3, cs, sn);
  M2d_make_translation(m4, a[0], a[1]);

  M2d_make_identity(m);
  M2d_mat_mult(m, m1, m);
  M2d_mat_mult(m, m2, m);
  M2d_mat_mult(m, m3, m);
  M2d_mat_mult(m, m4, m);
  M2d_mat_mult_points(rx,ry, m, rx,ry, 8);

						   

  G_rgb(0,1,1) ;
  G_fill_polygon(rx,ry,8) ;
  G_wait_key() ;
}
