//$Header$
//------------------------------------------------------------------------------
//                           PreliminaryOrbitDeterminationTest
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of test function of PreliminaryOrbitDeterminationTest that
* that are necesary for the proyect.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef PRELIMINARY_ORBIT_DETERMINATION_TEST
#define PRELIMINARY_ORBIT_DETERMINATION_TEST

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "doubler.h"
#include "angl.h"
#include "newtonnu.h"
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"
#include "Frac.h"
#include "timediff.h"
#include "Mjday.h"
#include "Position.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "gibbs.h"

void testUnit();
void testDoubler();
void testAngl();
void testNewtonnu();
void testR_x();
void testR_y();
void testR_z();
void testFrac();
void testTimediff();
void testMjday();
void testPosition();
void testMeanObliquity();
void testNutAngles();
void testGibbs();

#endif /* PRELIMINARY_ORBIT_DETERMINATION_TEST_H */