//$Header$
//------------------------------------------------------------------------------
//                                   PoleMatrix.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of PoleMatrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef POLEMATRIX_H
#define POLEMATRIX_H

#include "MatLabUtilites.h"
#include "R_y.h"
#include "R_x.h"

void PoleMatrix(double xp, double yp, double PoleMat[ROWS][COLS]);

#endif /* POLEMATRIX_H */