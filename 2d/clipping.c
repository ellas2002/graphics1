//Sutherland-Hodgman Polygon Clipping
#include "../FPToolkit.c"
#include <stdio.h>
#define MAXPOINTS 100

// saves the clicks the user inputs  
int click_and_save(double x[], double y[], int max_points){
  G_rgb(0,1,0);
  G_fill_rectangle(0,0,20,10);
  
  double v[2];
  int n = 0;
  while(1){
      G_wait_click(v);
      
      if(v[0] < 20 && v[1] < 10) {break;}
      x[n] = v[0]; y[n] = v[1]; n++;
      G_fill_circle(v[0],v[1],2);
      
      if (n >= max_points) break;
    }
  
  return n;
}

// Checks if a point is inside the clipping edge
int inside(double x, double y, double x1, double y1, double x2, double y2) {
    return (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) >= 0;
}

// Compute intersection point of line through (x1,y1)(x2,y2) 
// and line through (x3,y3)(x4,y4)
void intersection(double x1, double y1, double x2, double y2, 
                 double x3, double y3, double x4, double y4,
                 double *x, double *y) {
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    *x = ((x1*y2 - y1*x2) * (x3 - x4) - (x1 - x2) * (x3*y4 - y3*x4)) / denom;
    *y = ((x1*y2 - y1*x2) * (y3 - y4) - (y1 - y2) * (x3*y4 - y3*x4)) / denom;
}

// Sutherland-Hodgman polygon clipping
void suthHodge(double polyX[], double polyY[], int polysize, 
               double clipX[], double clipY[], int clipsize,
               double *outputX, double *outputY, int *outputsize) {
    double inputX[MAXPOINTS], inputY[MAXPOINTS];
    int inputSize = polysize;
    
    // Copy initial polygon
    for(int i = 0; i < polysize; i++) {
        inputX[i] = polyX[i];
        inputY[i] = polyY[i];
    }
    
    // Clip against each edge of the clipping polygon
    for(int j = 0; j < clipsize; j++) {
        double clipEdgeX1 = clipX[j];
        double clipEdgeY1 = clipY[j];
        double clipEdgeX2 = clipX[(j+1) % clipsize];
        double clipEdgeY2 = clipY[(j+1) % clipsize];
        
        // Output points for this iteration
        double outputIterX[MAXPOINTS], outputIterY[MAXPOINTS];
        int outputIterSize = 0;
        
        // Check each point of the input polygon against this clipping edge
        for(int i = 0; i < inputSize; i++) {
            double currX = inputX[i];
            double currY = inputY[i];
            double nextX = inputX[(i+1) % inputSize];
            double nextY = inputY[(i+1) % inputSize];
            
            // First point
            if(inside(currX, currY, clipEdgeX1, clipEdgeY1, clipEdgeX2, clipEdgeY2)) {
                // If inside, add to output
                outputIterX[outputIterSize] = currX;
                outputIterY[outputIterSize] = currY;
                outputIterSize++;
                
                // If line crosses, add intersection
                if(!inside(nextX, nextY, clipEdgeX1, clipEdgeY1, clipEdgeX2, clipEdgeY2)) {
                    double intX, intY;
                    intersection(currX, currY, nextX, nextY, 
                                 clipEdgeX1, clipEdgeY1, clipEdgeX2, clipEdgeY2, 
                                 &intX, &intY);
                    outputIterX[outputIterSize] = intX;
                    outputIterY[outputIterSize] = intY;
                    outputIterSize++;
                }
            } else {
                // If first point is outside, check if line crosses
                if(inside(nextX, nextY, clipEdgeX1, clipEdgeY1, clipEdgeX2, clipEdgeY2)) {
                    double intX, intY;
                    intersection(currX, currY, nextX, nextY, 
                                 clipEdgeX1, clipEdgeY1, clipEdgeX2, clipEdgeY2, 
                                 &intX, &intY);
                    outputIterX[outputIterSize] = intX;
                    outputIterY[outputIterSize] = intY;
                    outputIterSize++;
                }
            }
        }
        
        // Update input for next iteration
        inputSize = outputIterSize;
        for(int i = 0; i < inputSize; i++) {
            inputX[i] = outputIterX[i];
            inputY[i] = outputIterY[i];
        }
    }
    
    // Return final clipped polygon
    *outputsize = inputSize;
    for(int i = 0; i < inputSize; i++) {
        outputX[i] = inputX[i];
        outputY[i] = inputY[i];
    }
}

int main()
{
  G_init_graphics(800, 800);
  G_rgb(0.3, 0.3, 0.3);
  G_clear();
  
  double ax[MAXPOINTS], ay[MAXPOINTS];  // Subject polygon 
  double cx[MAXPOINTS], cy[MAXPOINTS];  // Clipping polygon
  
  G_rgb(1,0,0);
  int polysize = click_and_save(ax, ay, MAXPOINTS);
  G_rgb(0,1,0.5);
  G_fill_polygon(ax, ay, polysize);
  
  G_rgb(1,1,0);
  int clipsize = click_and_save(cx, cy, MAXPOINTS);
  G_rgb(0,1,1);
  G_fill_polygon(cx, cy, clipsize);
  
  // Clipped polygon storage
  double clippedX[MAXPOINTS], clippedY[MAXPOINTS];
  int clippedsize = 0;
  
  // Perform polygon clipping
  suthHodge(ax, ay, polysize, cx, cy, clipsize, 
            clippedX, clippedY, &clippedsize);
  
  // Draw clipped polygon
  G_rgb(0,0,1);  // Blue for clipped polygon
  G_fill_polygon(clippedX, clippedY, clippedsize);
  
  int key;   
  key = G_wait_key();
  
  return 0;
}
