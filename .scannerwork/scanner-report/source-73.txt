//$Header$
//------------------------------------------------------------------------------
//                                   GHAMatrix
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of GHAMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "GHAMatrix.h"

//------------------------------------------------------------------------------
//  void GHAMatrix(double Mjd_UT1, double ** eop, double filas,double columnas,
//				   double GHAmat[ROWS][COLS]);
//------------------------------------------------------------------------------
/**
* Transformation from true equator and equinox to Earth equator
*            and Greenwich meridian system
*
* @param  -
* Input:
*   Mjd_UT1   Modified Julian Date UT1
* @return -
* Output:
*   GHAmat    Greenwich Hour Angle matrix
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void GHAMatrix(double Mjd_UT1, double ** eop, double filas,double columnas,double GHAmat[ROWS][COLS]) {
	double gast_result = gast(Mjd_UT1,eop,filas,columnas);
    R_z(gast_result,GHAmat);
}