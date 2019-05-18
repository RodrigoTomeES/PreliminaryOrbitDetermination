//$Header$
//------------------------------------------------------------------------------
//                                   hgibbs
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of hgibbs function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "hgibbs.h"

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

void hgibbs(double r1[], double r2[], double r3[], double MJD1, double MJD2, double MJD3, double v2[], double *theta, double *theta1, double *copa, char error[]) {
	strcpy(error,"ok");
	*theta = 0;
	*theta1= 0;
	double magr1 = norm(r1);
	double magr2 = norm(r2);
	double magr3 = norm(r3);

	for (int i = 0; i < 3; i++) {
	    v2[i] = 0;
	}

	double tolangle= 0.01745329251994;
	double dt21= (MJD2-MJD1)*86400;
	double dt31= (MJD3-MJD1)*86400;
	double dt32= (MJD3-MJD2)*86400;

	double p[3], pn[3], r1n[3];
	crossVector(r2,r3,p);
	unit(p,pn);
	unit(r1,r1n);
	*copa = asin(dot(pn,r1n));

	if ( fabs(dot(r1n,pn)) > 0.017452406 ) {
		strcpy(error,"not coplanar");
	}

	*theta  = angl(r1,r2);
	*theta1 = angl(r2,r3);

	if ( (*theta > tolangle) || (*theta1 > tolangle) ) {
	    strcpy(error,"angl > 1ø");
	}

	double term1= -dt32*( 1/(dt21*dt31) + GM_Earth/(12*magr1*magr1*magr1) );
	double term2= (dt32-dt21)*( 1/(dt21*dt32) + GM_Earth/(12*magr2*magr2*magr2) );
	double term3=  dt21*( 1/(dt32*dt31) + GM_Earth/(12*magr3*magr3*magr3) );

	double r1_aux[3], r2_aux[3], r3_aux[3], sum_aux[3];
	// v2 =  term1*r1 + term2* r2 + term3* r3;
	multiplicacionVectorPorEscalar(r1,term1,r1_aux);
	multiplicacionVectorPorEscalar(r2,term2,r2_aux);
	multiplicacionVectorPorEscalar(r3,term3,r3_aux);
	sumaVectores(r1_aux, r2_aux, sum_aux);
	sumaVectores(r3_aux, sum_aux, v2);
}