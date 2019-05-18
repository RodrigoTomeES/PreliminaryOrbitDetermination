#include "PoleMatrix.h"

void PoleMatrix(double xp, double yp, double PoleMat[ROWS][COLS]){
    double Ry[ROWS][COLS];
    double Rx[ROWS][COLS];

    R_y(-xp,Ry);
    R_x(-yp,Rx);
    crossMatrix(Ry,Rx,PoleMat);
}