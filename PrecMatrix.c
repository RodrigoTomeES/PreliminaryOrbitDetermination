//$Header$
//------------------------------------------------------------------------------
//                                   PrecMatrix
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of PrecMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "PrecMatrix.h"

//------------------------------------------------------------------------------
//  void PrecMatrix(double Mjd_1, double Mjd_2, double PrecMat[ROWS][COLS])
//------------------------------------------------------------------------------
/**
* Precession transformation of equatorial coordinates
*
* @param  -
* Inputs:
*   Mjd_1     Epoch given (Modified Julian Date TT)
*   MjD_2     Epoch to precess to (Modified Julian Date TT)
* @return -
* Output:
*   PrecMat   Precession transformation matrix
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void PrecMatrix(double Mjd_1, double Mjd_2, double PrecMat[ROWS][COLS]) {
    double T = (Mjd_1-MJD_J2000)/36525;
    double dT= (Mjd_2-Mjd_1)/36525;

    double zeta = ( (2306.2181+(1.39656-0.000139*T)*T)+((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
    double z = zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
    double theta = ( (2004.3109-(0.85330+0.000217*T)*T)-((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;

    double Rz [ROWS][COLS];
    double Ry [ROWS][COLS];
    double Rz2 [ROWS][COLS];

    R_z(-z,Rz);
    R_y(theta,Ry);
    R_z(-zeta, Rz2);

    double aux[ROWS][COLS];
    crossMatrix(Rz,Ry,aux);
    crossMatrix(aux,Rz2,PrecMat);
}