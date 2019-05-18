//$Header$
//------------------------------------------------------------------------------
//                           anglesdr.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of anglesdr function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef ANGLESDR_H
#define ANGLESDR_H

#include "MatLabUtilites.h"
#include "doubler.h"
#include "lambert_gooding.h"
#include "stdio.h"

void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1, double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3, double * rsite1, double * rsite2, double * rsite3,double **  r2,double **  v2);

#endif /* ANGLESDR_H */