#ifndef IERS_H
#define IERS_H

#include "MatLabUtilites.h"

void IERS(double ** eop, double filas,double columnas,double Mjd_UTC, char interp, double * UT1_UTC, double * TAI_UTC, double * x_pole, double * y_pole, double * ddpsi, double * ddeps);

#endif /* IERS_H */