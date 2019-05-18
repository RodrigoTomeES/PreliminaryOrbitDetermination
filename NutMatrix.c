//$Header$
//------------------------------------------------------------------------------
//                                   NutMatrix
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of NutMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "NutMatrix.h"
#define TAM 3

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
void NutMatrix(double Mjd_TT, double NutMat[TAM][TAM]) {
    // Mean obliquity of the ecliptic
    double ep = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles(Mjd_TT, &dpsi, &deps);

    // Transformation from mean to true equator and equinox
    // NutMat = R_x(-ep-deps)*R_z(-dpsi)*R_x(+ep);
    double r_x[TAM][TAM], r_z[TAM][TAM], r_x2[TAM][TAM];

    R_x(-ep-deps,r_x);
    muestraMatriz(r_x);
    R_z(-dpsi,r_z);
    muestraMatriz(r_z);
    R_x(+ep,r_x2);
    muestraMatriz(r_x2);

    double aux[TAM][TAM];
    crossMatrix(r_x, r_z, aux);
    crossMatrix(aux, r_x2, NutMat);
}