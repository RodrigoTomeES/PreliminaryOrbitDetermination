//$Header$
//------------------------------------------------------------------------------
//                           gibbs.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of gibbs function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef GIBBS_H
#define GIBBS_H

#include <string.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "angl.h"

//------------------------------------------------------------------------------
//  void gibbs(double r1[], double r2[], double r3[], double res_vector[])
//------------------------------------------------------------------------------
/**
* Performs the gibbs method of orbit determination. this method
* determines the velocity at the middle point of the 3 given
* position vectors.
*
* @param - double r1[], double r2[], double r3[], double res_vector[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void gibbs(double r1[], double r2[], double r3[], double v2[], double *theta, double *theta1, double *copa, char error[]);

#endif /* GIBBS_H */