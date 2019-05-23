//$Header$
//------------------------------------------------------------------------------
//                                   Angl
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic implementation of Matlab functions that are not in C.
*
* @note
*/
//------------------------------------------------------------------------------

#include "angl.h"

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
double angl(double vector1[], double vector2[]) {
    double small     = 0.00000001;
    double undefined = 999999.1;
    double magv1 = norm(vector1);
    double magv2 = norm(vector2);
    double temp, theta;

    if (magv1*magv2 > (small*small)) {
        temp= dot(vector1,vector2)/(magv1*magv2);
        if (fabs(temp) > 1) {
            temp= sign(temp);
        }
        theta= acos(temp);
    } else {
        theta= undefined;
    }

    return theta;
}