#ifndef GAST_H
#define GAST_H

#include "MatLabUtilites.h"
#include "IERS.h"
#include "timediff.h"
#include "gmst.h"
#include "EqnEquinox.h"

double gast(double Mjd_UT1, double ** eop, double filas,double columnas);

#endif /* GAST_H */