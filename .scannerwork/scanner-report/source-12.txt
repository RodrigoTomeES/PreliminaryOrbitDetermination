//$Header$
//------------------------------------------------------------------------------
//                                   EqnEquinox
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of EqnEquinox function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "NutAngles.h"
#include "MeanObliquity.h"
#include "math.h"
#include "EqnEquinox.h"

//------------------------------------------------------------------------------
//  double EqnEquinox(double x);
//------------------------------------------------------------------------------
/**
* Computation of the equation of the equinoxes
*
* @param  -
* Input:
*   Mjd_TT    Modified Julian Date (Terrestrial Time)
* @return -
* Output:
*    EqE      Equation of the equinoxes
* @exception - none
* @see - none
* @note - The equation of the equinoxes dpsi*cos(eps) is the right ascension of
*   	  the mean equinox referred to the true equator and equinox and is equal
*   	  to the difference between apparent and mean sidereal time.
*/
//------------------------------------------------------------------------------
double EqnEquinox(double x){
    double dpsi,deps;
    NutAngles(x,&dpsi,&deps);

    double res = dpsi*cos(MeanObliquity(x));
    return res;
}