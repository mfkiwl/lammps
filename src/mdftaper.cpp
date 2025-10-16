#include <math.h>

/* MDF taper function between rmin and rmax as in GULP */

void mdftaper(double r, double rmin, double rmax, double &f, double &df)
{
  double x, dx, t1, t2;

  if (r <= rmin) {
    f = 1.0;
    df = 0.0;
  } else if (r >= rmax) {
    f = 0.0;
    df = 0.0;
  } else {
    x = (r-rmin)/(rmax-rmin);
    dx = 1.0/(rmax-rmin);
    t1 = 1.0 - x;
    t2 = 1.0 + 3.0*x + 6.0*pow(x,2);
    f = pow(t1,3)*t2;
    df = (-3.0*pow(t1,2)*t2 + pow(t1,3)*(3.0+12.0*x))*dx;
  }
}
