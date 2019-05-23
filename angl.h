//$Header$
//------------------------------------------------------------------------------
//                           angl.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of Matlab functions that are not in C.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef ANGL_H
#define ANGL_H

#include "MatLabUtilites.h"

//------------------------------------------------------------------------------
//  double norm(double vector[])
//------------------------------------------------------------------------------
/**
* Execute the test for MatLabUtilites.c
* This function show the results of the tests.
*
* @param  - double vector[]
* @return - double
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double angl(double vector1[], double vector2[]);

#endif /* ANGL_H */