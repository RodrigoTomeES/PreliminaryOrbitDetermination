#ifndef LAMBERTGOODING_H
#define LAMBERTGOODING_H

#include <stdlib.h>
#include <stdbool.h>
#include "MatLabUtilites.h"

void tlamb(double m,double q,double qsqfm1,double x,double n,double * t,double * dt,double * d2t,double * d3t);

double d8rt(double x);

#endif /* LAMBERTGOODING_H */