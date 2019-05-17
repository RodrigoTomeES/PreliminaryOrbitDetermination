#include "GHAMatrix.h"

void GHAMatrix(double Mjd_UT1, double ** eop, double filas,double columnas,double GHAmat[ROWS][COLS]){
	double gast_result = gast(Mjd_UT1,eop,filas,columnas);
	printf("%f\n", gast_result);
    R_z(gast_result,GHAmat);
}