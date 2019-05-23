//$Header$
//------------------------------------------------------------------------------
//                                   PoleMatrix
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of PoleMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "PoleMatrix.h"

//------------------------------------------------------------------------------
//  void PoleMatrix(double xp, double yp, double PoleMat[ROWS][COLS]
//------------------------------------------------------------------------------
/**
* Transformation from pseudo Earth-fixed to Earth-fixed
*             coordinates for a given date
*
* @param  -
* Input:
*   Pole coordinte(xp,yp)
* @return -
* Output:
*   PoleMat   Pole matrix
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void PoleMatrix(double xp, double yp, double PoleMat[ROWS][COLS]) {
    double Ry[ROWS][COLS];
    double Rx[ROWS][COLS];

    R_y(-xp,Ry);
    R_x(-yp,Rx);
    crossMatrix(Ry,Rx,PoleMat);
}