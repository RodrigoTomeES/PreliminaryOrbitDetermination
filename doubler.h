//$Header$
//------------------------------------------------------------------------------
//                                   doubler.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of doubler function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef DOUBLER_H
#define DOUBLER_H

#include "MatLabUtilites.h"

//------------------------------------------------------------------------------
//  void doubler(double cc1, double cc2, double magrsite1, double magrsite2,
//              double magrlin, double magr2in, double los1[], double los2[],
//              double los3[], double rsite1[],double rsite2[], double rsite3[],
//              double t1, double t3, char direct, double r2 [], double r3 [],
//              double * f1,double *f2,double *q1,double * magr1, double * magr2,
//              double *a , double * deltae32)
//------------------------------------------------------------------------------
/**
* This rountine accomplishes the iteration work for the double-r angles
*
* @param  - double cc1, double cc2, double magrsite1, double magrsite2, double magrlin, double magr2in, double los1[], double los2[], double los3[], double rsite1[],double rsite2[], double rsite3[],double t1, double t3, char direct, double r2 [], double r3 [], double * f1,double *f2,double *q1,double * magr1, double * magr2,double *a , double * deltae32
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magrlin, double magr2in, double los1[], double los2[], double los3[], double rsite1[],double rsite2[], double rsite3[], double t1, double t3, char direct, double r2 [], double r3 [], double * f1,double *f2,double *q1,double * magr1, double * magr2,double *a , double * deltae32);

#endif /* DOUBLER_H */