#include "GHAMatrix.h"

void GHAMatrix(double Mjd_UT1, double ** eop, double filas,double columnas,double GHAmat[ROWS][COLS]){
    R_z(gast(Mjd_UT1,eop,filas,columnas),GHAmat);
}