//$Header$
//------------------------------------------------------------------------------
//                           NutMatrix.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of NutMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef NUT_MATRIX_H
#define NUT_MATRIX_H

#include "MatLabUtilites.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

void NutMatrix(double Mjd_TT, double NutMat[3][3]);

#endif /* NUT_MATRIX_H */