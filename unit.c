//$Header$
//------------------------------------------------------------------------------
//                                   Unit
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

#include "unit.h"

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
void unit(double vector[], double unit_vector[]) {
    double small = 0.000001;
    double magv = norm(vector);

    if ( magv > small ) {
        for (int i = 0; i < 3; i++) {
            unit_vector[i] = vector[i]/magv;
        }
    } else {
        for (int i = 0; i < 3; i++) {
            unit_vector[i] = 0.0;
        }
    }
}