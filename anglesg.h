//$Header$
//------------------------------------------------------------------------------
//                           anglesg.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of anglesg function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef ANGLESG_H
#define ANGLESG_H

#include "MatLabUtilites.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "lambert_gooding.h"
#include "rv2coe.h"

//------------------------------------------------------------------------------
//  void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1,
//               double Delta2, double Delta3, double JD1, double JD2, double JD3,
//               double RS1[], double RS2[], double RS3[], double R2[], double V2[])
//------------------------------------------------------------------------------
/**
* Solves the problem of orbit determination using three
*               optical sightings.
*
* @param -
*  Inputs:         description               range/units
*    rtasc1      - right ascension #1            rad
*    rtasc2      - right ascension #2            rad
*    rtasc3      - right ascension #3            rad
*    decl1       - declination #1                rad
*    decl2       - declination #2                rad
*    decl3       - declination #3                rad
*    jd1         - julian date of 1st sighting   days from 4713 bc
*    jd2         - julian date of 2nd sighting   days from 4713 bc
*    jd3         - julian date of 3rd sighting   days from 4713 bc
*    RS1         - ijk site1 position vector     [m]
*    RS2         - ijk site2 position vector     [m]
*    RS3         - ijk site3 position vector     [m]
* @return -
*  Outputs:
*    r            - ijk position vector at t2     m
*    v            - ijk velocity vector at t2     m/s
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1, double Delta2, double Delta3, double JD1, double JD2, double JD3, double RS1[], double RS2[], double RS3[], double R2[], double V2[]);

#endif /* ANGLESG_H */