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
#define PI 3.141592653589793

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
void newtonnu(double ecc, double nu, double result[]) {
    double e0= 999999.9;
    double m = 999999.9;
    double small = 0.00000001;

    //--------------------------- circular ------------------------
    if ( fabs(ecc) < small  ){
        m = nu;
        e0= nu;
    } else {
        //---------------------- elliptical -----------------------
        if ( ecc < 1.0-small  ){
            double sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
            double cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
            e0 = atan2( sine,cose );
            m = e0 - ecc*sin(e0);

        } else {
            //-------------------- hyperbolic  --------------------
            if ( ecc > 1.0 + small){
                if((ecc > 1.0)&(fabs(nu) + 0.00001 < PI-acos(1.0 /ecc))){
                    double sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                    e0  = asinh( sine );
                    m   = ecc*sinh(e0) - e0;
                }
            } else{
                // ----------------- parabolic ---------------------
                if ( fabs(nu) < 168.0*PI/180.0 ){

                    e0= tan( nu*0.5);
                    m = e0 + (e0*e0*e0)/3.0;
                }
            }
        }
    }

    if(ecc < 1.0){
        m = drem(m,PI*2.0);
        if (m < 0.0){
            m = m + (2.0 *PI);
        }
        e0 = drem(e0, 2.0 *PI);
    }

    result[0] = e0;
    result[1] = m;
}
