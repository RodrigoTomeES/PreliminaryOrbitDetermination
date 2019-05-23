//$Header$
//------------------------------------------------------------------------------
//                                   lambert_gooding.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of lambert_gooding function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef LAMBERT_GOODING_H
#define LAMBERT_GOODING_H

#include <stdlib.h>
#include <stdbool.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "stdio.h"

//------------------------------------------------------------------------------
//  void tlamb(double m, double q, double qsqfm1, double x, double n, double *t,
//             double *dt, double *d2t, double *d3t)
//------------------------------------------------------------------------------
/**
* Lambert Gooding use this method
*
* @param  - double m, double q, double qsqfm1, double x, double n
* @return - double *t, double *dt, double *d2t, double *d3t
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void tlamb(double m,double q,double qsqfm1,double x,double n,double * t,double * dt,double * d2t,double * d3t);

//------------------------------------------------------------------------------
//  double d8rt(double x)
//------------------------------------------------------------------------------
/**
* Lambert Gooding use this method
*
* @param  - double x
* @return - double
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double d8rt(double x);

//------------------------------------------------------------------------------
//  void xlamb(double m,double q,double qsqfm1,double tin, double * n,double * x,
//             double * xpl)
//------------------------------------------------------------------------------
/**
* Lambert Gooding use this method
*
* @param  - double m,double q,double qsqfm1,double tin
* @return - double * n,double * x,double * xpl
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void xlamb(double m,double q,double qsqfm1,double tin, double * n,double * x,double * xpl);

//------------------------------------------------------------------------------
//  void vlamb(double gm,double * r1,double * r2,double th,double tdelt,
//             double * n,double ** vri,double ** vti,double ** vrf,double ** vtf);
//------------------------------------------------------------------------------
/**
* Lambert Gooding use this method
*
* @param  - double gm, double * r1,double * r2,double th,double tdelt
* @return - double * n,double ** vri,double ** vti,double ** vrf,double ** vtf
* @exception - none
* @see - none
* @note -
*    https://stackoverflow.com/questions/33132384/c-function-returning-pointer-and-a-dynamic-array
*    You can modify your code to set a pointer passed into your function by address, like this:
*
*    void create(int n, int** res) {
*        *res = malloc(n*sizeof(int));
*    }
*
*    Here is how you call this function now:
*
*    int *nn
*    int n = 5;
*    create(n, &nn);
*/
//------------------------------------------------------------------------------
void vlamb(double gm,double * r1,double * r2,double th,double tdelt,double * n,double ** vri,double ** vti,double ** vrf,double ** vtf);

//------------------------------------------------------------------------------
//  void lambert_gooding(double * r1,double * r2,double tof,double mu,
//                      double long_way,double multi_revs,double **  v1,
//                      double **  v2)
//------------------------------------------------------------------------------
/**
* Lambert's problem using Gooding's method
*
* @param  -
*  r1            first cartesian position [km]
*  r2            second cartesian position [km]
*  tof           time of flight [sec]
*  mu            gravity parameter [km^3/s^2]
*  long_way      when true, do "long way" (>pi) transfers
*  multi_revs    maximum number of multi-rev solutions to compute
* @return -
*  v1            vector containing 3d arrays with the cartesian components
*                of the velocities at r1
*  v2            vector containing 3d arrays with the cartesian components
*                of the velocities at r0
* @exception - none
* @see - none
* @note -
*  References:
*  1. R. H, Gooding. "[A procedure for the solution of Lambert's orbital
*     boundary-value problem](http://adsabs.harvard.edu/abs/1990CeMDA..48..145G)"
*     Celestial Mechanics and Dynamical Astronomy,
*     vol. 48, no. 2, 1990, p. 145-165.
*  2. A. Klumpp, "Performance Comparision of Lambert and Kepler Algorithms",
*     JPL Interoffice Memorandum, 314.1-0426-ARK, Jan 2, 1991.
*     [Zip](http://derastrodynamics.com/docs/lambert_papers_v1.zip)
*/
//------------------------------------------------------------------------------
void lambert_gooding(double * r1,double * r2,double tof,double mu,double long_way,double multi_revs,double **  v1,double **  v2);

#endif /* LAMBERT_GOODING_H */