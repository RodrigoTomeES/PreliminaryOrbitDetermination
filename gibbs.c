//$Header$
//------------------------------------------------------------------------------
//                                   gibbs
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of gibbs function.
*
* @note
*/
//------------------------------------------------------------------------------

#include <string.h>
#include "gibbs.h"
#include "unit.h"
#include "angl.h"
#define TAM 3

//------------------------------------------------------------------------------
//  void gibbs(double r1[], double r2[], double r3[], double res_vector[])
//------------------------------------------------------------------------------
/**
* Performs the gibbs method of orbit determination. this method
* determines the velocity at the middle point of the 3 given
* position vectors.
*
* @param - double r1[], double r2[], double r3[], double res_vector[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------

void gibbs(double r1[], double r2[], double r3[], double v2[], double *theta, double *theta1, double *copa, char error[]) {
	double small = 0.00000001;
	*theta = 0;
	*theta1 = 0;
	zeros_vector(v2);
	strcpy(error, "ok");

	double magr1 = norm(r1);
	double magr2 = norm(r2);
	double magr3 = norm(r3);

	double p[TAM],q[TAM],w[TAM];
	crossVector(r2,r3,p);
	crossVector(r3,r1,q);
	crossVector(r1,r2,w);

	double pn[TAM], r1n[TAM];
	unit(p,pn);
	unit(r1,r1n);
	*copa = asin(dot(pn,r1n));

	if ( fabs(dot(r1n,pn)) > 0.017452406 ) {
		strcpy(error, "not coplanar");
	}

	double d[TAM];
	double aux_d[TAM];
	// d = p + q + w;
	sumaVectores(p,q,aux_d);
	sumaVectores(aux_d,w,d);

	double magd = norm(d);

	double n[TAM];
	double magr1_aux[TAM], magr2_aux[TAM], magr3_aux[TAM], n_aux[TAM];
	// n = magr1*p + magr2*q + magr3*w;
	multiplicacionVectorPorEscalar(p,magr1,magr1_aux);
	multiplicacionVectorPorEscalar(q,magr2,magr2_aux);
	multiplicacionVectorPorEscalar(w,magr3,magr3_aux);
	sumaVectores(magr1_aux, magr2_aux, n_aux);
	sumaVectores(n_aux, magr3_aux, n);

	double magn = norm(n);
	double nn[TAM], dn[TAM];
	unit(n,nn);
	unit(d,dn);

	// -------------------------------------------------------------
	// determine if  the orbit is possible. both d and n must be in
	// the same direction, and non-zero.
	// -------------------------------------------------------------
	if ( ( fabs(magd)<small ) || ( fabs(magn)<small ) || ( dot(nn,dn) < small ) ){
		strcpy(error, "impossible");
	} else{
		*theta  = angl(r1,r2);
		*theta1 = angl(r2,r3);

		// ----------- perform gibbs method to find v2 -----------
/*		double r1mr2[TAM], r3mr1[TAM], r2mr3[TAM];
		restaVectores(magr1,magr2,r1mr2);
		restaVectores(magr3,magr1,r3mr1);
		restaVectores(magr2,magr3,r2mr3);*/
		double r1mr2 = magr1-magr2;
		double r3mr1 = magr3-magr1;
		double r2mr3 = magr2-magr3;

		double s[TAM];
		double r1mr2_aux[TAM], r3mr1_aux[TAM], r2mr3_aux[TAM], s_aux[TAM];
		// s  = r1mr2*r3 + r3mr1*r2 + r2mr3*r1;
		multiplicacionVectorPorEscalar(r3,r1mr2,r1mr2_aux);
		multiplicacionVectorPorEscalar(r2,r3mr1,r3mr1_aux);
		multiplicacionVectorPorEscalar(r1,r2mr3,r2mr3_aux);
		sumaVectores(r1mr2_aux, r3mr1_aux, s_aux);
		sumaVectores(s_aux, r2mr3_aux, s);

		double b[TAM];
		crossVector(d,r2,b);
		double l  = sqrt(GM_Earth/(magd*magn));
		double tover2 = l/magr2;

		double b2_aux[TAM], s2_aux[TAM];
		// v2 = tover2 * b + l * s;
		multiplicacionVectorPorEscalar(b,tover2,b2_aux);
		multiplicacionVectorPorEscalar(s,l,s2_aux);
		sumaVectores(b2_aux, s2_aux, v2);
	}
}