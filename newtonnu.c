//$Header$
//------------------------------------------------------------------------------
//                                   Newtonnu
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of newtonnu function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "newtonnu.h"

//------------------------------------------------------------------------------
//  double norm(double vector[])
//------------------------------------------------------------------------------
/**
*  newtonnu: Solves keplers equation when the true anomaly is known.
*    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
*    the parabolic limit at 168 is arbitrary. the hyperbolic anomaly is also
*    limited. the hyperbolic sine is used because it's not double valued.
*
* @param -
*  Inputs:         description                    range / units
*    ecc         - eccentricity                   0.0  to
*    nu          - true anomaly                   -2pi to 2pi rad
* @return -
*  Outputs:
*    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 deg
*    m           - mean anomaly                   0.0  to 2pi rad       151.7425 deg
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void newtonnu(double ecc, double nu, double * e0, double * m) {
    double e0V= 999999.9;
    double mV = 999999.9;
    double small = 0.00000001;

    //--------------------------- circular ------------------------
    if ( fabs(ecc) < small  ){
        mV = nu;
        e0V= nu;
    } else {
        //---------------------- elliptical -----------------------
        if ( ecc < 1.0-small  ){
            double sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
            double cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
            e0V = atan2( sine,cose );
            mV = e0V - ecc*sin(e0V);

        } else {
            //-------------------- hyperbolic  --------------------
            if ( ecc > 1.0 + small){
                if((ecc > 1.0)&(fabs(nu) + 0.00001 < M_PI-acos(1.0 /ecc))){
                    double sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                    e0V  = asinh( sine );
                    mV   = ecc*sinh(e0V) - e0V;
                }
            } else{
                // ----------------- parabolic ---------------------
                if ( fabs(nu) < 168.0*M_PI/180.0 ){

                    e0V= tan( nu*0.5);
                    mV = e0V + (e0V*e0V*e0V)/3.0;
                }
            }
        }
    }

    if(ecc < 1.0){
        mV = drem(mV,M_PI*2.0);
        if (mV < 0.0){
            mV = mV + (2.0 *M_PI);
        }
        e0V = drem(e0V, 2.0 *M_PI);
    }

    *e0 = e0V;
    *m = mV;
}
