//$Header$
//------------------------------------------------------------------------------
//                                   MeanObliquity
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of MeanObliquity function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "MeanObliquity.h"

//------------------------------------------------------------------------------
//  double MeanObliquity(double Mjd_TT)
//------------------------------------------------------------------------------
/**
* Computes the mean obliquity of the ecliptic
*
* @param - Mjd_TT - Modified Julian Date (Terrestrial Time)
* @return - MOblq - Mean obliquity of the ecliptic
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT-MJD_J2000)/36525;

    double MOblq = Rad*(23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600);

    return MOblq;
}