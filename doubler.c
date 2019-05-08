#include "doubler.h"
#define TAM 3

void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magrlin, double magr2in, double los1[], double los2[], double los3[], double rsite1[],double rsite2[], double rsite3[], double t1, double t3, char direct, double * magr1, double * magr2, double r3 []){
    double mu = 398600.4418e9;
    double rho1 = (-cc1 + sqrt(cc1*cc1-4*(magrsite1*magrsite1-magrlin*magrlin)))/2.0;
    double rho2 = (-cc2 + sqrt(cc2*cc2-4*(magrsite2*magrsite2-magr2in*magr2in)))/2.0;

    double r1[TAM];
    double r2[TAM];


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


}