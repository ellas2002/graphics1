#include "FPToolkit.c"

int main() {
  int numpoints, i = 0, n = 0 ;
  double x[4500], y[4500] ;
  int numpolys ;
  int psize[9000] ;
  int con[9000][31] ;
  double red[9000], green[9000], blue[9000] ;

  scanf("%d",&numpoints) ;
  for (i=0 ; i < numpoints ; i++) {
    scanf("%lf %lf", &x[i],&y[i]) ;
  }
  scanf("%d", &numpolys) ;
  for (i=0;i<numpolys;i++) {
    scanf("%d", &psize[i]) ;
    for (n=0;n<psize[i];n++) {
      scanf("%d", &con[i][n]) ;
    }
  }
  for (i=0;i<numpolys;i++) {
    scanf("%lf %lf %lf", &red[i], &green[i], &blue[i]) ;
  }


  G_init_graphics(800,800) ;

  double a[9000], b[9000] ;
  
  for(i=0;i<numpolys;i++) {
    for(n=0;n<psize[i];n++) {
      a[n] = x[con[i][n]] ;
      b[n] = y[con[i][n]] ;
    }
    G_rgb(red[i],green[i],blue[i]) ;
    G_fill_polygon(a,b,psize[i]) ;
  }

  int key ;
  key = G_wait_key() ;

}
