//$Header$
//------------------------------------------------------------------------------
//                                   doubler
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a implementation of doubler function.
*
* @note
*/
//------------------------------------------------------------------------------

#include "doubler.h"
#include <stdio.h>
#define TAM 3

//------------------------------------------------------------------------------
//  void doubler(double cc1, double cc2, double magrsite1, double magrsite2,
//              double magrlin, double magr2in, double los1[], double los2[],
//              double los3[], double rsite1[],double rsite2[], double rsite3[],
//              double t1, double t3, char direct, double r2 [], double r3 [],
//              double * f1,double *f2,double *q1,double * magr1, double * magr2,
//              double *a , double * deltae32)
//------------------------------------------------------------------------------
/**
* This rountine accomplishes the iteration work for the double-r angles
*
* @param  - double vector[]
* @return - double
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magrlin, double magr2in, double los1[], double los2[], double los3[], double rsite1[],double rsite2[], double rsite3[],double t1, double t3, char direct, double r2 [], double r3 [], double * f1,double *f2,double *q1,double * magr1, double * magr2,double *a , double * deltae32) {

    double mu = 398600.4418e9;
    double rho1 = (-cc1 + sqrt(cc1*cc1-4*(magrsite1*magrsite1-magrlin*magrlin)))/2.0;
    double rho2 = (-cc2 + sqrt(cc2*cc2-4*(magrsite2*magrsite2-magr2in*magr2in)))/2.0;

    double r1[TAM];


    double aux[TAM];
    double aux1[TAM];

    multiplicacionVectorPorEscalar(los1,rho1,aux);
    sumaVectores(aux,rsite1,r1);

    multiplicacionVectorPorEscalar(los2,rho2,aux);
    sumaVectores(aux,rsite2,r2);

    *magr1=norm(r1);
    *magr2=norm(r2);

    double w [TAM];

    if(direct=='y'){
        crossVector(r1,r2,aux);
        devisionVectorPorEscalar(aux,(*magr1)*(*magr2),w);
    }else{
        crossVector(r1,r2,aux);
        devisionVectorPorEscalar(aux,(*magr1)*(*magr2),aux1);
        opuestoVector(aux1,w);
    }

    double rho3;
    rho3=-dot(rsite3,w)/dot(los3,w);

    multiplicacionVectorPorEscalar(los3,rho3,aux);
    sumaVectores(aux,rsite3,r3);

    double magr3 = norm(r3);

    double cosdv21 = dot(r2,r1)/((*magr2)*(*magr1));
    double sindv21 = sqrt(1-cosdv21*cosdv21);
    double dv21 = atan2(sindv21,cosdv21);

    double cosdv31 = dot(r3,r1)/((magr3)*(*magr1));
    double sindv31 = sqrt(1-cosdv31*cosdv31);
    double dv31 = atan2(sindv31,cosdv31);

    double cosdv32 = dot(r3,r2)/((magr3)*(*magr2));
    double sindv32 = sqrt(1-cosdv32*cosdv32);

    double pi = 3.14159265358979;

    double c1 ,c3,p;

    if(dv31>pi){
        c1 = ((*magr2)*sindv32)/((*magr1)*sindv31);
        c3= ((*magr2)*sindv21)/((magr3)*sindv31);
        p=(c1*(*magr1)+c3*magr3-(*magr2))/(c1+c3-1);
    }else{
        c1 = ((*magr1)*sindv31)/((*magr2)*sindv32);
        c3= ((*magr1)*sindv21)/((magr3)*sindv32);
        p=(c3*(magr3)-c1*(*magr2)+(*magr1))/(-c1+c3+1);
    }

    double ecosv1,ecosv2,ecosv3;
    ecosv1=p/(*magr1)-1;
    ecosv2=p/(*magr2)-1;
    ecosv3=p/(magr3)-1;

    double esinv2;

    if(dv21!=pi){
        esinv2=(-cosdv21*ecosv2+ecosv1)/sindv21;
    }else{
        esinv2 = (cosdv32*ecosv2-ecosv3)/sindv31;
    }

    double e;
    e = sqrt(ecosv2*ecosv2+esinv2*esinv2);
    *a = p/(1-e*e);

    double n,s,c;

    double sinde32,cosde32,sinde21,cosde21,deltae21,deltam32,deltam12,sindh32,sindh21,deltah32,deltah21;
    if(e*e<1){
        n=sqrt(mu/((*a)*(*a)*(*a)));
        s=(*magr2)/p*sqrt(1-e*e)*esinv2;
        c=(*magr2)/p*(e*e+ecosv2);

        sinde32 = magr3/sqrt((*a)*p)*sindv32-magr3/p*(1-cosdv32)*s;
        cosde32 = 1-(*magr2)*magr3/((*a)*p)*(1-cosdv32);
        *deltae32 = atan2(sinde32,cosde32);

        sinde21 = (*magr1)/sqrt((*a)*p)*sindv21+(*magr1)/p*(1-cosdv21)*s;
        cosde21 = 1-(*magr2)*(*magr1)/((*a)*p)*(1-cosdv21);
        deltae21 = atan2(sinde21,cosde21);

        deltam32 = (*deltae32)+2*s*(sin((*deltae32)/2))*(sin((*deltae32)/2))-c*sin((*deltae32));
        deltam12 = -deltae21+2*s*(sin(deltae21/2))*(sin(deltae21/2))+c*sin(deltae21);
    }else{
        n=sqrt(mu/-((*a)*(*a)*(*a)));

        s = (*magr2)/p*sqrt(e*e-1)*esinv2;
        c = (*magr2)/p*(e*e+ecosv2);

        sindh32 = magr3/sqrt(-(*a)*p)*sindv32-magr3/p*(1-cosdv32)*s;
        sindh21 = (*magr1)/sqrt(-(*a)*p)*sindv21+(*magr1)/p*(1-cosdv21)*s;

        deltah32 = log( sindh32 + sqrt(sindh32*sindh32 +1) );
        deltah21 = log( sindh21 + sqrt(sindh21*sindh21 +1) );

        deltam32 = -deltah32+2*s*(sinh(deltah32/2))*(sinh(deltah32/2))+c*sinh(deltah32);
        deltam12 = deltah21+2*s*(sinh(deltah21/2))*(sinh(deltah21/2))-c*sinh(deltah21);
       *deltae32 = deltah32;
    }

    *f1=t1-deltam12/n;
    *f2=t3-deltam32/n;

    *q1=sqrt((*f1)*(*f1)+(*f2)*(*f2));
}