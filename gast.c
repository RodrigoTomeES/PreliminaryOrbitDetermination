//$Header$
//------------------------------------------------------------------------------
//                                   Gast
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of Gast function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "gast.h"

//------------------------------------------------------------------------------
//  double gast(double Mjd_UT1, double ** eop, double filas,double columnas)
//------------------------------------------------------------------------------
/**
* Greenwich Apparent Sidereal Time
*
* @param  -
* Input:
*   Mjd_UT1   Modified Julian Date UT1
* @return -
* Output:
*   gstime    GAST in [rad]
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double gast(double Mjd_UT1, double ** eop, double filas,double columnas) {
    double UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps;

    IERS(eop,filas,columnas,Mjd_UT1,'l',&UT1_UTC,&TAI_UTC,&x_pole,&y_pole,&ddpsi,&ddeps);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_UTC,TAI_UTC,&UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    double Mjd_UTC = Mjd_UT1-UT1_UTC/86400;
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;

    double gstime = matlab_mod(gmst(Mjd_UT1)+EqnEquinox(Mjd_TT), pi2);
    return gstime;
}