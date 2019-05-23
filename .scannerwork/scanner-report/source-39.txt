//$Header$
//------------------------------------------------------------------------------
//                           hgibbs.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of hgibbs function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef HGIBBS_H
#define HGIBBS_H

#include <string.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "angl.h"

//------------------------------------------------------------------------------
//  void hgibbs(double r1[], double r2[], double r3[], double MJD1, double MJD2,
//				double MJD3, double v2[], double *theta, double *theta1,
//				double *copa, char error[]);
//------------------------------------------------------------------------------
/**
* implements the herrick-gibbs approximation for orbit
*          determination, and finds the middle velocity vector for the 3
*          given position vectors.
*
* @param -
*  inputs:
*    r1          - ijk position vector #1         m
*    r2          - ijk position vector #2         m
*    r3          - ijk position vector #3         m
*    MJD1        - julian date of 1st sighting    days from 4713 bc
*    MJD2        - julian date of 2nd sighting    days from 4713 bc
*    MJD3        - julian date of 3rd sighting    days from 4713 bc
* @return -
*  outputs:
*    v2          - ijk velocity vector for r2     m/s
*    theta       - angl between vectors           rad
*    error       - flag indicating success        'ok',...
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void hgibbs(double r1[], double r2[], double r3[], double MJD1, double MJD2, double MJD3, double v2[], double *theta, double *theta1, double *copa, char error[]);

#endif /* HGIBBS_H */