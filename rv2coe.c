//$Header$
//------------------------------------------------------------------------------
//                                   rv2coe
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of rv2coe function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "rv2coe.h"

//------------------------------------------------------------------------------
//  void rv2coe(double r[], double v[], double *p, double *a, double *ecc,
//				double *incl, double *omega, double *argp, double *nu, double *m,
//				double *arglat, double *truelon, double *lonper)
//------------------------------------------------------------------------------
/**
* Finds the classical orbital elements given the geocentric
*          equatorial position and velocity vectors.
*
* @param -
*	Inputs:         description                    range / units
*	    r           - ijk position vector            m
*	    v           - ijk velocity vector            m/s
* @return -
*  Outputs:
*    p           - semilatus rectum               m
*    a           - semimajor axis                 m
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------

void rv2coe(double r[], double v[], double *p, double *a, double *ecc, double *incl, double *omega, double *argp, double *nu, double *m, double *arglat, double *truelon, double *lonper) {
double mu = 398600.4418e9;
double small = 1e-10;
double undefined = 999999.1;
double magr= norm(r);
double magv= norm(v);
// ------------------  find h n and e vectors   ----------------
double hbar[3];
crossVector( r,v,hbar );
double magh= norm( hbar );
if ( magh > small )
	double nbar[3];
    nbar[0]= -hbar[1];
    nbar[1]=  hbar[0];
    nbar[2]=   0.0;
    double magn = norm( nbar );
    double c1 = magv*magv - mu /magr;
    double rdotv= dot( r,v );
    double ebar[3];
    for (int i = 0; i < 3; i++) {
        ebar[i]= (c1*r[i] - rdotv*v[i])/mu;
    }

    *ecc = norm( ebar );

    // ------------  find a e and semi-latus rectum   ----------
    double sme= ( magv*magv*0.5  ) - ( mu /magr );

    if ( fabs( sme ) > small ) {
        *a= -mu  / (2.0 *sme);
    } else {
        *a= infinite; //WTF
    }
    *p = magh*magh/mu;

    // -----------------  find inclination   -------------------
    double hk= hbar[2]/magh;
    *incl= acos( hk );

    // --------  determine type of orbit for later use  --------
    // ------ elliptical, parabolic, hyperbolic inclined -------
    char typeorbit[]= 'ei';
    if ( *ecc < small ){
        // ----------------  circular equatorial ---------------
        if  (*incl<small) || (fabs(*incl-pi)<small){
            strcpy(typeorbit,"ce");
          }else{
            // --------------  circular inclined ---------------
            strcpy(typeorbit,"ci");
        }
      }else{
        // - elliptical, parabolic, hyperbolic equatorial --
        if  (*incl<small) || (fabs(*incl-pi)<small){
            strcpy(typeorbit,"ee");
        }
    }

    // ----------  find longitude of ascending node ------------
    if ( magn > small )
        temp= nbar(1) / magn;
        if ( fabs(temp) > 1.0  )
            temp= sign(temp);
        end
        *omega= acos( temp );
        if ( nbar(2) < 0.0  )
            *omega= 2*pi - *omega;
        end
      else
        *omega= undefined;
    end

    // ---------------- find argument of perigee ---------------
    if ( strcmp(typeorbit,'ei') )
        *argp = angl( nbar,ebar);
        if ( ebar(3) < 0.0  )
            *argp= 2*pi - *argp;
        end
      else
        *argp= undefined;
    end

    // ------------  find true anomaly at epoch    -------------
    if ( typeorbit(1:1) == 'e' )
        *nu =  angl( ebar,r);
        if ( rdotv < 0.0  )
            *nu= 2*pi - *nu;
        end
      else
        *nu= undefined;
    end

    // ----  find argument of latitude - circular inclined -----
    if ( strcmp(typeorbit,'ci') )
        *arglat = angl( nbar,r );
        if ( r(3) < 0.0  )
            *arglat= 2*pi - *arglat;
        end
        *m = *arglat;
      else
        *arglat= undefined;
    end

    // -- find longitude of perigee - elliptical equatorial ----
    if  ( *ecc>small ) && (strcmp(typeorbit,'ee'))
        temp= ebar(1)/(*ecc);
        if ( fabs(temp) > 1.0  )
            temp= sign(temp);
        end
        *lonper= acos( temp );
        if ( ebar(2) < 0.0  )
            *lonper= 2*pi - *lonper;
        end
        if ( *incl > halfpi )
            *lonper= 2*pi - *lonper;
        end
      else
        *lonper= undefined;
    end

    // -------- find true longitude - circular equatorial ------
    if  ( magr>small ) && ( strcmp(typeorbit,'ce') )
        temp= r(1)/magr;
        if ( fabs(temp) > 1.0  )
            temp= sign(temp);
        end
        *truelon= acos( temp );
        if ( r(2) < 0.0  )
            *truelon= 2*pi - *truelon;
        end
        if ( *incl > halfpi )
            *truelon= 2*pi - *truelon;
        end
        m = *truelon;
      else
        *truelon= undefined;
    end

    // ------------ find mean anomaly for all orbits -----------
    if ( typeorbit(1:1) == 'e' ){
        [e,m] = newtonnu(*ecc,*nu );
    }
else{
    *p    = undefined;
    a    = undefined;
    *ecc  = undefined;
    *incl = undefined;
    *omega= undefined;
    *argp = undefined;
    *nu   = undefined;
    *m    = undefined;
    *arglat = undefined;
    *truelon= undefined;
    *lonper = undefined;
}
}