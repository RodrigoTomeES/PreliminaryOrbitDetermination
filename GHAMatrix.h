//$Header$
//------------------------------------------------------------------------------
//                                   GHAMAtrix.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of GHAMAtrix function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef GHAMATRIX_H
#define GHAMATRIX_H

#include "MatLabUtilites.h"
#include "R_z.h"
#include "gast.h"

void GHAMatrix(double Mjd_UT1, double ** eop, double filas,double columnas,double GHAmat[ROWS][COLS]);

#endif /* GHAMATRIX_H */