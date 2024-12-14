#include "FPToolkit.c"
#include "M2d_matrix_tools.c"

double jx[5] = {0, 40, 35, 10,  0} ;
double jy[5] = {0,  0,  8,  8, 15} ;
double ojx[5] = {0, 40, 35, 10,  0} ;
double ojy[5] = {0,  0,  8,  8, 15} ;


int main()
{
  double m[3][3], m1[3][3], m2[3][3], m3[3][3], m4[3][3], m5[3][3];
  double dx, dy, length, a;
  int unit = 0;
  double x = 0;double y = 500;
  int key, frame = 0;
  char frame_str[100];

  
  G_init_graphics(600,600) ;

  dx = 500 - 0;
  dy = 0 - 500;
  length = sqrt((dx*dx) + (dy*dy));
  a = atan2(dy, dx) ;
  
  double flip_state=0;
  
  M2d_make_scaling(m1, 1, 1);
  M2d_make_rotation(m2, a);
  M2d_make_translation(m4, x, y);

  
  M2d_make_identity(m);
  M2d_mat_mult(m, m1, m2);
  M2d_mat_mult(m, m4, m);

  M2d_mat_mult_points(jx, jy, m, jx, jy, 5); 
    
  G_rgb(0,0,0) ;  G_clear() ;
  G_rgb(0,0,1) ;  G_fill_polygon(jx,jy,5) ;
  G_rgb(1,1,0) ;  G_draw_string("any key to continue", 250,10) ;
  G_rgb(1,0,1) ;  G_line(0,500,500,0) ;
  
  while(1){
     key = G_wait_key();
        
     if (key == 'q' || key == 'Q') {
        printf("Quitting... Total frames: %d\n", frame);
        break;
     }
	
     frame++;
        
     if (frame % 50 == 0) {
            flip_state = !flip_state;
     }

     x += 2.0 * (dx/length);
     y += 2.0 * (dy/length);
        
     M2d_make_identity(m);

     if (flip_state) {
       M2d_make_translation(m3, x, y);
       M2d_make_scaling(m5, 1, -1);
       M2d_make_rotation(m2, M_PI/4);
       M2d_mat_mult(m, m5, m2);      
     }

     else {
       M2d_make_translation(m3, x, y);
       M2d_make_scaling(m1, 1, 1);
       M2d_make_rotation(m2, a);
       M2d_mat_mult(m, m1, m2);
     }

     if(frame % 250 == 0){
       M2d_make_translation(m3, y, x);
       M2d_make_scaling(m1, 1, 1);
       M2d_make_rotation(m2, a);
       M2d_mat_mult(m, m1, m2);
     }

     M2d_mat_mult(m, m3, m);
        

     M2d_mat_mult_points(jx, jy, m, ojx, ojy, 5);
        

     G_rgb(0,0,0); G_clear();
     G_rgb(0,0,1); G_fill_polygon(jx, jy, 5);
     G_rgb(1,1,0);
     sprintf(frame_str, "Frame: %d - Press any key to continue, 'q' to quit\n", frame);
     printf("Frame: %d - Press any key to continue, 'q' to quit\n", frame);
     G_draw_string(frame_str, 250,10);
     G_rgb(1,0,1); G_line(0,500,500,0);
    }
  
  G_close() ;
  return 0;

}
