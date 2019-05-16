#include "IERS.h"

void IERS(double ** eop, double filas,double columnas,double Mjd_UTC, char interp, double * UT1_UTC, double * TAI_UTC, double * x_pole, double * y_pole, double * ddpsi, double * ddeps){
    double pi = 3.14159265358979;
    
    //double Arcs = 3600*180/pi;
    
    if(interp=='1'){
        double mj = floor(Mjd_UTC);
        double nop = columnas;

        double * preeop;
        double * nexteop;
        for(int i=0;i<nop;i++){
            if(mj== eop[4][i]){
                for(int j=0;j<filas;j++){
                    preeop[j]=eop[j][i];
                    nexteop[j]=eop[j][i+1];
                }
                break;
            }
        }

        double mfme = 1440*(Mjd_UTC-mj);
        double fixf=mfme/1440;

        *UT1_UTC = preeop[6]+(nexteop[6]-preeop[6])*fixf;
        *TAI_UTC=preeop[12];

        *x_pole=preeop[4]+(nexteop[4]-preeop[4])*fixf;
        *y_pole=preeop[5]+(nexteop[5]-preeop[5])*fixf;

        *ddpsi=preeop[8]+(nexteop[8]-preeop[8])*fixf;
        *ddeps=preeop[9]+(nexteop[9]-preeop[9])*fixf;

        *x_pole=((*x_pole)/Arcs);
        *y_pole=((*y_pole)/Arcs);

        *ddpsi=((*ddpsi)/Arcs);
        *ddeps=((*ddeps)/Arcs);
    }else if(interp == 'n'){
        double mj = floor(Mjd_UTC);
        double nop = columnas;

        double * preeop;
        for(int i=0;i<nop;i++){
            if(mj== eop[4][i]){
                for(int j=0;j<filas;j++){
                    preeop[j]=eop[j][i];
                }
                break;
            }
        }

        *UT1_UTC = preeop[6];
        *TAI_UTC=preeop[12];

        *x_pole=preeop[4]/Arcs;
        *y_pole=preeop[5]/Arcs;

        *ddpsi=preeop[8]/Arcs;
        *ddeps=preeop[9]/Arcs;

    }
}