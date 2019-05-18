//$Header$
//------------------------------------------------------------------------------
//                                   anglesg
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of anglesg function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "anglesg.h"

//------------------------------------------------------------------------------
//  void anglesg(double r1[], double r2[], double r3[], double res_vector[])
//------------------------------------------------------------------------------
/**
* Solves the problem of orbit determination using three
*               optical sightings.
*
* @param -
*  Inputs:         description               range/units
*    rtasc1      - right ascension #1            rad
*    rtasc2      - right ascension #2            rad
*    rtasc3      - right ascension #3            rad
*    decl1       - declination #1                rad
*    decl2       - declination #2                rad
*    decl3       - declination #3                rad
*    jd1         - julian date of 1st sighting   days from 4713 bc
*    jd2         - julian date of 2nd sighting   days from 4713 bc
*    jd3         - julian date of 3rd sighting   days from 4713 bc
*    RS1         - ijk site1 position vector     [m]
*    RS2         - ijk site2 position vector     [m]
*    RS3         - ijk site3 position vector     [m]
* @return -
*  Outputs:
*    r            - ijk position vector at t2     m
*    v            - ijk velocity vector at t2     m/s
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------

void gibbs(double r1[], double r2[], double r3[], double v2[], double *theta, double *theta1, double *copa, char error[]) {
Mu = 398600.4418e9;
Rad = 180/pi;

R1 = zeros(3,1);
R2 = zeros(3,1);
R3 = zeros(3,1);

//---------- set middle to 0, find deltas to others -----------
Tau1 = (JD1-JD2)*86400;
Tau3 = (JD3-JD2)*86400;

//----------------  Find Line of SIGHT vectors  ---------------
L1(1) = cos(Delta1)*cos(Alpha1);
L1(2) = cos(Delta1)*sin(Alpha1);
L1(3) = sin(Delta1);

L2(1) = cos(Delta2)*cos(Alpha2);
L2(2) = cos(Delta2)*sin(Alpha2);
L2(3) = sin(Delta2);

L3(1) = cos(Delta3)*cos(Alpha3);
L3(2) = cos(Delta3)*sin(Alpha3);
L3(3) = sin(Delta3);

//------------- Find L matrix and determinant -----------------
//--------- Called LMatI since it is only used for determ -----
for i = 1:3
    LMatIi(i,1) =L1(i);
    LMatIi(i,2) =L2(i);
    LMatIi(i,3) =L3(i);
    RSMat(i,1) =RS1(i);
    RSMat(i,2) =RS2(i);
    RSMat(i,3) =RS3(i);
end

D = det(LMatIi);

//------------------ Now assign the inverse -------------------
LMatI(1,1) = ( L2(2)*L3(3)-L2(3)*L3(2))/D;
LMatI(2,1) = (-L1(2)*L3(3)+L1(3)*L3(2))/D;
LMatI(3,1) = ( L1(2)*L2(3)-L1(3)*L2(2))/D;
LMatI(1,2) = (-L2(1)*L3(3)+L2(3)*L3(1))/D;
LMatI(2,2) = ( L1(1)*L3(3)-L1(3)*L3(1))/D;
LMatI(3,2) = (-L1(1)*L2(3)+L1(3)*L2(1))/D;
LMatI(1,3) = ( L2(1)*L3(2)-L2(2)*L3(1))/D;
LMatI(2,3) = (-L1(1)*L3(2)+L1(2)*L3(1))/D;
LMatI(3,3) = ( L1(1)*L2(2)-L1(2)*L2(1))/D;

LIR = LMatI*RSMat;

//------------ Find f and g series at 1st and 3rd obs ---------
// *      speed by assuming circ sat vel for uDot here ??
// *      some similartities in 1/6t3t1 ...
//--- keep separated this time ----
a1  = Tau3/(Tau3 - Tau1);
a1u = (Tau3*((Tau3-Tau1)*(Tau3-Tau1) - Tau3*Tau3 ))/(6.0*(Tau3 - Tau1));
a3  = -Tau1 / (Tau3 - Tau1);
a3u = -(Tau1*((Tau3-Tau1)*(Tau3-Tau1) - Tau1*Tau1 ))/(6.0*(Tau3 - Tau1));

//--- Form initial guess of r2 ----
D1 = LIR(2,1)*a1 - LIR(2,2) + LIR(2,3)*a3;
D2 = LIR(2,1)*a1u + LIR(2,3)*a3u;

//------- Solve eighth order poly NOT same as LAPLACE ---------
L2DotRS= dot(L2,RS2);
magRS2 = norm(RS2);
Poly( 1)=  1.0;  // r2^8th variable!!!!!!!!!!!!!!
Poly( 2)=  0.0;
Poly( 3)=  -(D1*D1 + 2.0*D1*L2DotRS + magRS2^2);
Poly( 4)=  0.0;
Poly( 5)=  0.0;
Poly( 6)=  -2.0*Mu*(L2DotRS*D2 + D1*D2);
Poly( 7)=  0.0;
Poly( 8)=  0.0;
Poly( 9)=  -Mu*Mu*D2*D2;
Poly(10)=  0.0;
Poly(11)=  0.0;
Poly(12)=  0.0;
Poly(13)=  0.0;
Poly(14)=  0.0;
Poly(15)=  0.0;
Poly(16)=  0.0;
rootarr = roots(Poly);

//------------------ Select the correct root ------------------
BigR2 = 0.0;

for j= 1:15
    if ( rootarr(j) > BigR2 ) && ( isreal(rootarr(j)) )
        BigR2 = rootarr(j);
    end
end

//------------ Solve matrix with u2 better known --------------
u = Mu/(BigR2*BigR2*BigR2);

c1 = a1+a1u*u;
c3 = a3+a3u*u;
CMat(1,1) = -c1;
CMat(2,1) = 1.0;
CMat(3,1) = -c3;
RhoMat = LIR*CMat;

// Rhoold1 =  RhoMat(1,1)/c1;
Rhoold2 = -RhoMat(2,1);
// Rhoold3 =  RhoMat(3,1)/c3;

//-------- Loop through the refining process ------------  for WHILE () DO
// for ll = 1:3
Rho2 = 999999e3;
ll = 0;
while ( (abs(Rhoold2-Rho2)>1e-12) && (ll<=2) )
    ll = ll+1;
    Rho2 = Rhoold2;  // reset now that inside while loop
    //---------- Now form the three position vectors ----------
    for i = 1:3
        R1(i)=  RhoMat(1,1)*L1(i)/c1 + RS1(i);
        R2(i)= -RhoMat(2,1)*L2(i)    + RS2(i);
        R3(i)=  RhoMat(3,1)*L3(i)/c3 + RS3(i);
    end

    [~,theta,theta1,copa,error] = gibbs(R1,R2,R3);
    if ( ~strcmp(error,'ok') && (copa < 1/Rad) )
        //--- HGibbs to get middle vector ----
        [~,theta,theta1,copa,error] = hgibbs(R1,R2,R3,JD1,JD2,JD3);
    end
    [~,V2] = lambert_gooding(R1',R2',(JD2-JD1)*86400,Mu,0,1);
    [p,a,ecc,incl,omega,argp,Nu,m,u,l,ArgPer] = rv2coe(R2,V2);
    magR2 = norm(R2);

    if ( ll <= 2 )
        //--- Now get an improved estimate of the f and g series --
        //       .or. can the analytic functions be found now??
        U = Mu/(magR2^3);
        RDot = dot(R2,V2)/magR2;
        UDot = (-3.0*Mu*RDot)/(magR2^4);

        TauSqr= Tau1*Tau1;
        f1 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau1 ...
            + (1.0/24.0) * U*U*TauSqr*TauSqr ...
            + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau1;
        g1 = Tau1 - (1.0/6.0)*U*Tau1*TauSqr - (1.0/12.0) * ...
            UDot*TauSqr*TauSqr ...
            + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau1 ...
            + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
        TauSqr = Tau3*Tau3;
        f3 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau3 ...
            + (1.0/24.0) * U*U*TauSqr*TauSqr ...
            + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau3;
        g3 = Tau3 - (1.0/6.0)*U*Tau3*TauSqr - (1.0/12.0) * ...
            UDot*TauSqr*TauSqr ...
            + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau3 ...
            + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
    else
        //-------- Now use exact method to find f and g -----------
        Theta = angl(R1,R2);
        Theta1 = angl(R2,R3);
        magR1 = norm(R1);
        magR3 = norm(R3);

        f1 = 1.0 - ( (magR1*(1.0 - cos(Theta))/p ) );
        g1 = ( magR1*magR2*sin(-theta) ) / sqrt(p);  // - ANGLE because backwards!!
        f3 = 1.0 - ( (magR3*(1.0 - cos(Theta1))/p ) );
        g3 = ( magR3*magR2*sin(theta1) )/sqrt(p);
    end

    c1 =  g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);

    //----- Solve for all three ranges via matrix equation ----
    CMat(1,1) = -c1;
    CMat(2,1) = 1.0;
    CMat(3,1) = -c3;
    RhoMat = LIR*CMat;

//     Rhoold1 =  RhoMat(1,1)/c1;
    Rhoold2 = -RhoMat(2,1);
//     Rhoold3 =  RhoMat(3,1)/c3;
    //----------------- Check for convergence -----------------
end   // DO WHILE the ranges are still changing

//---------------- Find all three vectors ri ------------------
for i= 1:3
    R1(i) =  RhoMat(1,1)*L1(i)/c1 + RS1(i);
    R2(i) = -RhoMat(2,1)*L2(i)    + RS2(i);
    R3(i) =  RhoMat(3,1)*L3(i)/c3 + RS3(i);
end
}