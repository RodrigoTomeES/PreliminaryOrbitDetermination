//$Header$ 
//------------------------------------------------------------------------------ 
//                           PreliminaryOrbitDeterminationTest
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination. 
// 
// Legal: MIT  License
// 
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27 
// 
/**  
* Provides a basic test for all function of Preliminary Orbit Determination Proyect.  
*  
* @note     
*/ 
//------------------------------------------------------------------------------

// Test MatLabUtilites
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "MatLabUtilites.h"
#include "unit.h"
#define EPSILON pow(10, -7)
// Epsilon es de 10⁻7 ya que es la precisión que deja ver matlab para el test

//Colores para los mensajes
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

double v1[] = {0.0,0.0,0.0};
double v2[]={4,-1,2};
double v3[]={2,-2,-1};

//------------------------------------------------------------------------------ 
//  int main()
//------------------------------------------------------------------------------ 
/**  
* Execute the test for Preliminary Orbit Determination Proyect
* This function show the results of the tests.
*  
* @param  - none
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/ 
//------------------------------------------------------------------------------
int main () {
    // Test unit
    printf("---- Test UNIT ----\n");

    double v1_unitario[3];
    double v2_unitario[3];
    double v3_unitario[3];

    unit(v1, v1_unitario);
    unit(v2, v1_unitario);
    unit(v3, v3_unitario);

    //Vector v1
    assert(fabs(v1_unitario[0] - 0.0) < EPSILON);
    assert(fabs(v1_unitario[1] - 0.0) < EPSILON);
    assert(fabs(v1_unitario[2] - 0.0) < EPSILON);


    //Vector v2
    assert(fabs(v2_unitario[0] - 0.8729) < EPSILON);
    assert(fabs(v2_unitario[1] - -0.2182) < EPSILON);
    assert(fabs(v2_unitario[2] - 0.4364) < EPSILON);

    //Vector v3
    assert(fabs(v3_unitario[0] - 0.6667) < EPSILON);
    assert(fabs(v3_unitario[1] - -0.6667) < EPSILON);
    assert(fabs(v3_unitario[2] - -0.3333) < EPSILON);

    printf(GREEN "---- Pass Test UNIT ----\n" RESET);
}