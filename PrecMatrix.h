//$Header$
//------------------------------------------------------------------------------
//                                   PrecMatrix.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of PrecMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef PRECMATRIX_H
#define PRECMATRIX_H

#include "MatLabUtilites.h"
#include "R_y.h"
#include "R_z.h"

void PrecMatrix(double Mjd_1, double Mjd_2, double PrecMat[ROWS][COLS]);

#endif /* PRECMATRIX_H */