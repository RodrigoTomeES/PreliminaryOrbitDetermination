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
    R_z(-dpsi,r_z);
    R_x(+ep,r_x2);

    double aux[TAM][TAM];
    multiplicacionMatrices(r_x, r_z, aux);
    multiplicacionMatrices(r_x2, aux, NutMat);
}

/*Mjd_TT =

          54977.6815585532


NutMat =

         0.999999997895152     -5.95287721488282e-05     -2.58073726784525e-05
      5.95281964998585e-05         0.999999997979423     -2.23057957973194e-05
      2.58087004629423e-05       2.2304259484005e-05         0.999999999418215*/