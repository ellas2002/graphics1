#include <stdio.h>
#include <math.h>

// STUDENTS NEED TO FILL IN SOME OF THESE

/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( 1 )          (1)

instead of (x',y',1) = (x,y,1) * M  

*/


//prints the matrix
int M2d_print_mat (double a[3][3])
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 




//be able to copy a matrix
int M2d_copy_mat (double a[3][3], double b[3][3])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           a[r][c] = b[r][c] ;
      }
  }

  return 1 ;
} 




//set up:  1 0 0
//         0 1 0
//         0 0 1
// for any matrix of size n x n, multiplyiny by I leaves A unchanged
int M2d_make_identity (double a[3][3])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
} 




//combined with scaling ad translating through matrix multiplication
int M2d_make_translation (double a[3][3], double dx, double dy)
{
  M2d_make_identity(a) ;
  a[0][2] =  dx ;  a[1][2] = dy ;  
  return 1 ;
}



//combined with scaling ad translating through matrix multiplication
int M2d_make_scaling (double a[3][3], double sx, double sy)
{
  // YOU NEED TO FILL THIS IN
  M2d_make_identity(a);
  a[0][0] = sx; a[1][1] = sy;
  return 1 ;
}




//combined with scaling ad translating through matrix multiplication
int M2d_make_rotation_cs (double a[3][3], double cs, double sn)
// this assumes cosine and sine are already known
{
  M2d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][1] = -sn ;
  a[1][0] =   sn ;  a[1][1] =  cs ;

  return 1 ;
}


//gotta add radians
int M2d_make_rotation (double a[3][3], double radians)
{
  double cs = cos(radians);
  double sn = sin(radians);
  
  M2d_make_rotation_cs(a, cs, sn) ;

  return 1;

}






//multipling matrices 
int M2d_mat_mult (double res[3][3], double a[3][3], double b[3][3])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// M2d_mat_mult(p,  p,q) or M2d_mat_mult(p,  q,p) or  M2d_mat_mult(p, p,p)
{
  // YOU NEED TO FILL THIS IN
  double temp[3][3] = {0};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      temp[i][j] = 0;
      for(int k = 0; k < 3; k++){
	temp[i][j] += a[i][k] * b[k][j];
      }
    }
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      res[i][j] = temp[i][j];
    }
  }
  
  return 1 ;
}



//function multipes a 2d point ! by a trasnformation matrix m and stores the resutls in the point P
int M2d_mat_mult_pt (double P[2],   double m[3][3], double Q[2])
// P = m*Q
// SAFE, user may make a call like M2d_mat_mult_pt (W, m,W) ;
{
  double u,v ;

  u = m[0][0]*Q[0] + m[0][1]*Q[1] + m[0][2] ;
  v = m[1][0]*Q[0] + m[1][1]*Q[1] + m[1][2] ;

  P[0] = u ;
  P[1] = v ;
  
  return 1 ;
}


// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like M2d_mat_mult_points (x,y, m, x,y, n) ;

//allows for two sepearte arrays
int M2d_mat_mult_points (double X[], double Y[],
                         double m[3][3],
                         double x[], double y[], int numpoints)
{
  double q[2], p[2];
  for(int i = 0; i <= numpoints - 1; i++){
    q[0] = x[i];
    q[1] = y[i];
    M2d_mat_mult_pt(p, m, q);

    X[i] = p[0];
    Y[i] = p[1];
   
  }
  
  return 1;
}







