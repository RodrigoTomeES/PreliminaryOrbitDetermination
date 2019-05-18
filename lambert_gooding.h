#ifndef LAMBERTGOODING_H
#define LAMBERTGOODING_H

#include <stdlib.h>
#include <stdbool.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "stdio.h"

void tlamb(double m,double q,double qsqfm1,double x,double n,double * t,double * dt,double * d2t,double * d3t);

double d8rt(double x);

void xlamb(double m,double q,double qsqfm1,double tin, double * n,double * x,double * xpl);
void vlamb(double gm,double * r1,double * r2,double th,double tdelt,double * n,double ** vri,double ** vti,double ** vrf,double ** vtf);

void lambert_gooding(double * r1,double * r2,double tof,double mu,double long_way,double multi_revs,double **  v1,double **  v2);

#endif /* LAMBERTGOODING_H */