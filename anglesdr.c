//$Header$
//------------------------------------------------------------------------------
//                                   anglesdr
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of anglesdr function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "angl.h"

//------------------------------------------------------------------------------
//  void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1,
//                double decl2, double decl3, double Mjd1, double Mjd2,
//                double Mjd3, double rsite1[], double rsite2[], double rsite3[],
//                double r2[], double v2[])
//------------------------------------------------------------------------------
/**
* Solves the problem of orbit determination using three
*
* @param  - double vector[]
* @return - double
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1, double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3, double rsite1[], double rsite2[], double rsite3[], double r2[], double v2[]) {
magr1in = 2.01*R_Earth;
magr2in = 2.11*R_Earth;
direct  = 'y';

tol    = 1e-8*R_Earth;
pctchg = 5e-6;

t1 = (Mjd1 - Mjd2)*86400;
t3 = (Mjd3 - Mjd2)*86400;

los1 = [cos(decl1)*cos(rtasc1) cos(decl1)*sin(rtasc1) sin(decl1)];
los2 = [cos(decl2)*cos(rtasc2) cos(decl2)*sin(rtasc2) sin(decl2)];
los3 = [cos(decl3)*cos(rtasc3) cos(decl3)*sin(rtasc3) sin(decl3)];

magr1old  = 99999e3;
magr2old  = 99999e3;
magrsite1 = norm(rsite1);
magrsite2 = norm(rsite2);
magrsite3 = norm(rsite3);

cc1 = 2*dot(los1,rsite1);
cc2 = 2*dot(los2,rsite2);

ll = 0;
while (abs(magr1in-magr1old) > tol && abs(magr2in-magr2old) > tol && ll<=3)
    ll = ll+1;
    [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler(cc1,cc2,magrsite1,...
     magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,...
     t3,direct);

//     [v2,v3] = lambert_gooding(r2',r3',(Mjd3-Mjd2)*86400,GM_Earth,0,1);
    f  = 1 - a/magr2*(1-cos(deltae32));
    g  = t3 - sqrt(a^3/GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - f*r2)/g;

    magr1o = magr1in;
    magr1in = (1+pctchg)*magr1in;
    deltar1 = pctchg*magr1in;
    [r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32] = doubler(cc1,cc2,...
     magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,...
     rsite3,t1,t3,direct);

    pf1pr1 = (f1delr1-f1)/deltar1;
    pf2pr1 = (f2delr1-f2)/deltar1;

    magr1in = magr1o;
    deltar1 = pctchg*magr1in;
    magr2o = magr2in;
    magr2in = (1+pctchg)*magr2in;
    deltar2 = pctchg*magr2in;
    [r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32] = doubler(cc1,cc2,...
     magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,...
     rsite3,t1,t3,direct);
    pf1pr2 = (f1delr2-f1)/deltar2;
    pf2pr2 = (f2delr2-f2)/deltar2;

    magr2in = magr2o;
    deltar2 = pctchg*magr2in;

    delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
    delta1 = pf2pr2*f1 - pf1pr2*f2;
    delta2 = pf1pr1*f2 - pf2pr1*f1;

    deltar1 = -delta1/delta;
    deltar2 = -delta2/delta;

    magr1old = magr1in;
    magr2old = magr2in;

    // may need to limit the amount of the correction
//     if abs(deltar1) > magr1in*pctchg
//         deltar1 = sign(deltar1)*magr1in*pctchg;
//     end
//     if abs(deltar2) > magr2in*pctchg
//         deltar2 = sign(deltar2)*magr2in*pctchg;
//     end

    magr1in = magr1in + deltar1;
    magr2in = magr2in + deltar2;
end

[r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler(cc1,cc2,magrsite1,...
 magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);

[v2,v3] = lambert_gooding(r2',r3',(Mjd3-Mjd2)*86400,GM_Earth,0,1);
f  = 1 - a/magr2*(1-cos(deltae32));
g  = t3 - sqrt(a^3/GM_Earth)*(deltae32-sin(deltae32));
v2 = (r3 - f*r2)/g;
}