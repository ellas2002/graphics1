#include "FPToolkit.c"
#include "M2d_matrix_tools.c"


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


int main()
{
   double a[100], b[100];  // Increased array size for safety
    double m1[3][3], m2[3][3], m3[3][3], m4[3][3], m[3][3];
    int nab;
    
    G_init_graphics(800, 800);
    G_rgb(0.3,0.3,0.3);
    G_clear();
    
    G_rgb(0,1,0);
    double p[2], q[2];
    
    // Create a polygon
    nab = click_and_save(a, b);
    G_polygon(a, b, nab);
    G_fill_polygon(a,b, nab);
    
    // Create reflection line
    G_wait_click(p);
    G_fill_circle(p[0], p[1], 2);
    G_wait_click(q);
    G_fill_circle(q[0], q[1], 2);
    
    G_rgb(0,1,0.5);
    G_line(p[0], p[1], q[0], q[1]);
    
    // Calculate reflection transformation
    double dx = q[0] - p[0];
    double dy = q[1] - p[1];
    double len = sqrt(dx*dx + dy*dy);
    double cs = dx/len;
    double sn = dy/len;
    
    // Build transformation matrix
    M2d_make_identity(m);
    M2d_make_translation(m1, -p[0], -p[1]);    // Move to origin
    M2d_mat_mult(m, m1, m);
    
    M2d_make_rotation_cs(m2, cs, -sn);         // Rotate to align with x-axis
    M2d_mat_mult(m, m2, m);
    
    M2d_make_scaling(m3, 1, -1);               // Reflect across x-axis
    M2d_mat_mult(m, m3, m);
    
    M2d_make_rotation_cs(m4, cs, sn);          // Rotate back
    M2d_mat_mult(m, m4, m);
    
    double m5[3][3];
    M2d_make_translation(m5, p[0], p[1]);      // Move back to original position
    M2d_mat_mult(m, m5, m);
    
    // Apply transformation to polygon points
    M2d_mat_mult_points(a, b, m, a, b, nab);
    
    // Draw reflected polygon
    G_rgb(0,1,0);
    G_polygon(a, b, nab);
    
    G_wait_key();
    return 0;
}
