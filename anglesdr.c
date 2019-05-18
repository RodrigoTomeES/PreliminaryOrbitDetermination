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

#include "anglesdr.h"

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
void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1, double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3, double * rsite1, double * rsite2, double * rsite3,double **  r2,double **  v2) {
    double magr1in = 2.01*R_Earth;
    double magr2in = 2.11*R_Earth;
    char direct  = 'y';

    double tol    = 1e-8*R_Earth;
    double pctchg = 5e-6;

    double t1 = (Mjd1 - Mjd2)*86400;
    double t3 = (Mjd3 - Mjd2)*86400;

    double los1[ROWS];
    los1[0]=cos(decl1)*cos(rtasc1);
    los1[1]=cos(decl1)*sin(rtasc1);
    los1[2]=sin(decl1);

    double los2[ROWS];
    los2[0]=cos(decl2)*cos(rtasc2);
    los2[1]=cos(decl2)*sin(rtasc2);
    los2[2]=sin(decl2);
     
    double los3[ROWS];
    los3[0]=cos(decl3)*cos(rtasc3);
    los3[1]=cos(decl3)*sin(rtasc3);
    los3[2]=sin(decl3);
     

    double magr1old  = 99999e3;
    double magr2old  = 99999e3;
    double magrsite1 = norm(rsite1);
    double magrsite2 = norm(rsite2);
    double magrsite3 = norm(rsite3);

    double cc1 = 2*dot(los1,rsite1);
    double cc2 = 2*dot(los2,rsite2);

    double aux1[ROWS];
    double aux2[ROWS];
    double aux3[ROWS];
    double v2v[ROWS];


    double r2v [ROWS]; double r3 [ROWS]; double f1, f2, q1, magr1, magr2, a, deltae32;
    double f1delr1, f2delr1, q2;
    double f1delr2,f2delr2, q3;
    double f,g;

    //printf("soy anglesdr\n");

    int ll = 0;
    while (fabs(magr1in-magr1old) > tol && fabs(magr2in-magr2old) > tol && ll<=3){
        //printf("bucleando\n");
        ll = ll+1;

        

        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2v,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);
        
        //     [v2v,v3] = lambert_gooding(r2v',r3',(Mjd3-Mjd2)*86400,GM_Earth,0,1);
        f  = 1 - a/magr2*(1-cos(deltae32));
        g  = t3 - sqrt((a*a*a)/GM_Earth)*(deltae32-sin(deltae32));
        


        //v2v = (r3 - f*r2v)/g;
        multiplicacionVectorPorEscalar(r2v,f,aux1);
        restaVectores(r3,aux1,aux2);
        devisionVectorPorEscalar(aux2,g,v2v);
        
        
        
        double magr1o = magr1in;
        magr1in = (1+pctchg)*magr1in;
        double deltar1 = pctchg*magr1in;

        
        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2v,r3,&f1delr1,&f2delr1,&q2,&magr1,&magr2,&a,&deltae32);
        
        double pf1pr1 = (f1delr1-f1)/deltar1;
        double pf2pr1 = (f2delr1-f2)/deltar1;
        
        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        double magr2o = magr2in;
        magr2in = (1+pctchg)*magr2in;
        double deltar2 = pctchg*magr2in;

        
        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2v,r3,&f1delr2,&f2delr2,&q3,&magr1,&magr2,&a,&deltae32);

        double pf1pr2 = (f1delr2-f1)/deltar2;
        double pf2pr2 = (f2delr2-f2)/deltar2;
        
        magr2in = magr2o;
        deltar2 = pctchg*magr2in;
        
        double delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        double delta1 = pf2pr2*f1 - pf1pr2*f2;
        double delta2 = pf1pr1*f2 - pf2pr1*f1;
        
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
    }

    doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,r2v,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);
    
    double * v2vaux;
    double * v3;
    
    lambert_gooding(r2v,r3,(Mjd3-Mjd2)*86400,GM_Earth,0,1,&v2vaux,&v3);
    
    f  = 1 - a/magr2*(1-cos(deltae32));
    g  = t3 - sqrt((a*a*a)/GM_Earth)*(deltae32-sin(deltae32));
    
    
    multiplicacionVectorPorEscalar(r2v,f,aux1);
    restaVectores(r3,aux1,aux2);
    devisionVectorPorEscalar(aux2,g,aux3);
    
    
    *r2 = (double *)calloc(3,sizeof(double));//v1 = zeros(3,n_solutions);
    *v2 = (double *)calloc(3,sizeof(double));//v2 = zeros(3,n_solutions);

    (*r2)[0] = r2v[0];
    (*r2)[1] = r2v[1];
    (*r2)[2] = r2v[2];
    
    (*v2)[0] = aux3[0];
    (*v2)[1] = aux3[1];
    (*v2)[2] = aux3[2];
    
}