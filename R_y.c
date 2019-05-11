//$Header$
//------------------------------------------------------------------------------
//                                   R_y
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of R_y function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "R_y.h"

//------------------------------------------------------------------------------
//  void R_y(double angle, double rotmat[][])
//------------------------------------------------------------------------------
/**
* Make the rot matrix
*
* @param - angle of rotation [rad]
* @return - rotation matrix
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void R_y(double angle, double rotmat[3][3]) {
    double C = cos(angle);
    double S = sin(angle);
    zeros(rotmat);

    rotmat[0][0] =   C;  rotmat[0][1] = 0.0;  rotmat[0][2] = -1.0*S;
    rotmat[1][0] = 0.0;  rotmat[1][1] = 1.0;  rotmat[1][2] =    0.0;
    rotmat[2][0] =   S;  rotmat[2][1] = 0.0;  rotmat[2][2] =      C;
}