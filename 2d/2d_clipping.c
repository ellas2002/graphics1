#include "../FPToolkit.c"

//Sutherland–Hodgman Algorithm

// saves the clicks the user inputs  
int click_and_save(double x[], double y[]){
  double v[2];
  int n = 0;

  while(1){
      G_wait_click(v); 
      if( v[1] < 40) {break;}

      x[n] = v[0];
      y[n] = v[1];
      n++;

      G_rgb(1, 1, 1);
      G_fill_circle(v[0],v[1],2);
    }
  
  return n;
 }



int in_out( double x[], double y[], int n,
	    double P[2], int I, int J)
{
  double mx = 0, my = 0;

  //get x and y coordiates of polygon
  for(int i = 0; i < n; i++){ 
    mx += x[i];
    my += y[i];
  }

  //calculate mean of x and y coordinates in polygon
  mx = mx / n;
  my = my / n;

  double a = (y[J] - y[I]); //change in y  between vertices
  double b = (x[J] - x[I]); // in x
  double c = (x[J] - x[I]) * y[I] - (y[J] - y[I]) * x[I]; //constant
  
  double signM = a * mx - b * my + c; //sign value for polygon centroid
  double signP = a * P[0] - b * P[1] + c; //for test point

  if (signM > 0 && signP < 0){return 0;} //on opposite sides of the line

  if (signM < 0 && signP > 0){return 0;} 


  return 1; //inside or on the line
}



//finds intersection of two lines
int intersect_2_lines (double A[2], double B[2],
		       double C[2], double D[2],
		       double intersection[2])
{
double x1 = A[0], y1 = A[1];
    double x2 = B[0], y2 = B[1];
    double x3 = C[0], y3 = C[1];
    double x4 = D[0], y4 = D[1];

    // find denominators
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    // check for parallel lines (denominator = 0)
    if (denom == 0) {
        return 0; // lines don't intersect
    }

    // find intersection point
    intersection[0] = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom;
    intersection[1] = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom;

    // check for intersection point in lines
    if (intersection[0] < fmin(x1, x2) || intersection[0] > fmax(x1, x2) ||
        intersection[0] < fmin(x3, x4) || intersection[0] > fmax(x3, x4) ||
        intersection[1] < fmin(y1, y2) || intersection[1] > fmax(y1, y2) ||
        intersection[1] < fmin(y3, y4) || intersection[1] > fmax(y3, y4)) {
        return 0; // intersection outside of lines
    }

    return 1; // lines intersect at point
}


//clips a polygon against a single edge of clipping polygon
int clip_line(double A[], double B[], int nab,
	      double x[], double y[], int nxy,
	      int a, int b)
{
  double P[2], intersection[2];
  double aa[2], bb[2], xx[2], yy[2];

  int inout1, inout2, i ,j;

  double tempx[1000], tempy[1000];
  int ntemp = 0;

  aa[0] = A[a]; aa[1] = B[a]; //sets up current clipping edge
  bb[0] = A[b]; bb[1] = B[b];


  G_rgb(0,1,1); //draws clipping polygon and clipping edge
  G_polygon(A, B, nab);
  G_rgb(1,0,1);
  G_line(A[a], B[a], A[b], B[b]);

  //iterates thorugh each edge of input polygon
  for( i = 0; i < nxy ; i++){
    j = (i + 1) % nxy;

    xx[0] = x[i]; xx[1] = y[i];
    yy[0] = x[j]; yy[1] = y[j];

    P[0] = x[i]; P[1] = y[i];
    inout1 = in_out(A, B, nab, P, a, b);

    P[0] = x[j]; P[1] = y[j];
    inout2 = in_out(A, B, nab, P, a, b); 

    //both points inside: keep second point
    if(inout1 == 1 && inout2 == 1){
      tempx[ntemp] = x[j];
      tempy[ntemp] = y[j];
      ntemp++;

      //first point inside, second outside: add intersecting point 
    }else if(inout1 == 1 && inout2 == 0){

      intersect_2_lines(aa, bb, xx, yy, intersection);
      G_rgb(0,1,0);

      G_fill_circle(intersection[0], intersection[1], 3);
      tempx[ntemp] = intersection[0];
      tempy[ntemp] = intersection[1];
      ntemp++;

      //first point outside, second inside: add intersection and second point
    }else if(inout1 == 0 && inout2 == 1){
      intersect_2_lines(aa, bb, xx, yy, intersection);
      G_rgb(0,1,0);

      G_fill_circle(intersection[0], intersection[1], 3);
      tempx[ntemp] = intersection[0];
      tempy[ntemp] = intersection[1];
      ntemp++;

      tempx[ntemp] = x[j];
      tempy[ntemp] = y[j];
      ntemp++;
    }
  }

  //draws original and clipped polygon
  G_rgb(1,0,0); G_polygon(x, y, nxy);
  G_rgb(0,1,1); G_polygon(A, B, nab);
  G_rgb(0,1,0); G_polygon(tempx, tempy, ntemp);

  //updates input polygon with clipped vertices
  for( i = 0; i < ntemp; i++){
    x[i] = tempx[i];
    y[i] = tempy[i];
  }

  return ntemp;
}




int main()
{
  double a[1000], b[1000];
  double x[1000], y[1000];
  double cx[1000], cy[1000];

  int nab, nxy, cxy;
  char q;

  G_init_graphics(800, 800);
  G_rgb(0,0,0);
  G_clear();
  G_rgb(1,0,0);
  G_fill_rectangle(0,0,800,40);

  //initial clipping polygon
  cxy = click_and_save(cx,cy);

  while(1){
    G_init_graphics(800, 800);
    G_rgb(0,0,0);
    G_clear();
    G_rgb(1,0,0);
    G_fill_rectangle(0,0,800,40);

    //draw polygons
    G_rgb(1,0,0);
    G_polygon(cx, cy, cxy);

    nab = click_and_save(a, b);
    G_rgb(0,1,1);
    G_polygon(a, b, nab);

    //nxy now has initial clipping points
    nxy = cxy;

    //copied to x and y
    for(int i = 0; i < cxy; i++){
      x[i] = cx[i];
      y[i] = cy[i];
    }

    //clip polygon
    for(int i = 0; i < nab; i++){
      int j = (i + 1) % nab;
      nxy = clip_line(a, b, nab, x, y, nxy, i, j);
    }

    //displays result
    G_init_graphics(800, 800);
    G_rgb(0,0,0);
    G_clear();
    G_rgb(1,0,0);
    G_fill_rectangle(0,0,800,40);


    G_rgb(0.3, 0, 0);
    G_polygon(cx, cy, cxy);
    G_rgb(0, 0.3, 0.3);
    G_polygon(a, b, nab);
    G_rgb(0, 1, 0);
    G_polygon(x, y, nxy);
   
    
    q = G_wait_key();
    if( q == 113){break;}  
  }

  return 0;
  
}
