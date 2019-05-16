//$Header$
//------------------------------------------------------------------------------
//                                   Frac
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of Frac function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "Frac.h"

//------------------------------------------------------------------------------
//  void Frac(double x, double res[][])
//------------------------------------------------------------------------------
/**
* Fractional part of a number (y=x-[x])
*
* @param - x
* @return - res
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double Frac(double x) {
	double y = x;
	printf("x %f \n",x);
	double res = floor(y);
	printf("floor %f \n",res);
	double caguedios =x-res;
	printf("caguen %f \n",caguedios);
	return caguedios;
	
}