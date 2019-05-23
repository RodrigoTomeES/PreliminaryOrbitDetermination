//$Header$
//------------------------------------------------------------------------------
//                                   Position
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of Position function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "Position.h"

//------------------------------------------------------------------------------
//  void Position(double lon, double lat, double h, double r[])
//------------------------------------------------------------------------------
/**
* Create a Position vector (r [m]) from geodetic coordinates
*
* @param - Longitude [rad], latitude [rad], altitude [m] and result vector
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void Position(double lon, double lat, double h, double r[]) {
    double R_equ = R_Earth;
    double f = f_Earth;

    double e2     = f*(2-f);  // Square of eccentricity
    double CosLat = cos(lat); // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ/sqrt(1-e2*SinLat*SinLat);

    r[0] = (       N+h)*CosLat*cos(lon);
    r[1] = (       N+h)*CosLat*sin(lon);
    r[2] = ((1-e2)*N+h)*SinLat;
}