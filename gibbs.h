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

#include "MatLabUtilites.h"

void gibbs(double r1[], double r2[], double r3[], double v2[], double *theta, double *theta1, double *copa, char error[]);

#endif /* GIBBS_H */