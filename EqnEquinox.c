#include "NutAngles.h"
#include "MeanObliquity.h"
#include "math.h"

double EqnEquinox(double x){
    double dpsi,deps;
    NutAngles(x,&dpsi,&deps);

    double res = dpsi*cos(MeanObliquity(x));
    return res;
}