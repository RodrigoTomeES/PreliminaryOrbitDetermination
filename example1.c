#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MatLabUtilites.h"
#include "Mjday.h"
#include "Position.h"
#include "timediff.h"
#include "IERS.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "anglesdr.h"

#define ROWS 3
#define COLS 3
#define TAM 3
//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
* Execute the test for Preliminary Orbit Determination Proyect
* This function show the results of the tests.
*
* @param  - none
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
int main () {
    // read Earth orientation parameters
    FILE* fid = fopen("eop19620101.txt","rt");

    int filas = 20026;
    int columnas = 13;
    int columna1, columna2, columna3, columna4, columna13;
    float columna5, columna6, columna7, columna8, columna9, columna10, columna11, columna12;

    if (fid == NULL){
        exit(EXIT_FAILURE);
    }

    double **eop;
    eop = (double **) malloc (filas*sizeof(double *));

    if (eop != NULL) {
        for (int i = 0; i < filas; i++) {
            eop[i] = (double *) malloc (columnas * sizeof(double));
            if (fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &columna1, &columna2, &columna3, &columna4, &columna5, &columna6, &columna7, &columna8, &columna9, &columna10, &columna11, &columna12, &columna13) == EOF) {
                break;
            }
            eop[i][0]  =  columna1;
            eop[i][1]  =  columna2;
            eop[i][2]  =  columna3;
            eop[i][3]  =  columna4;
            eop[i][4]  =  columna5;
            eop[i][5]  =  columna6;
            eop[i][6]  =  columna7;
            eop[i][7]  =  columna8;
            eop[i][8]  =  columna9;
            eop[i][9]  = columna10;
            eop[i][10] = columna11;
            eop[i][11] = columna12;
            eop[i][12] = columna13;

            //printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }
    } else {
        printf("Es null");
    }

    fclose(fid);

    // read observations
    fid = fopen("sat1.txt","rt");
    printf("\n");
    int filasOBS = 3;
    int columnasOBS = 3;
    int Y, M, D, h, m;
    float s, rtasc, decl;

    if (fid == NULL){
        exit(EXIT_FAILURE);
    }

    double **obs;
    obs = (double **) malloc (filasOBS*sizeof(double *));

    if (obs != NULL) {
        for (int i = 0; i < filasOBS; i++) {
            obs[i] = (double *) malloc (columnasOBS * sizeof(double));
            if (fscanf(fid,"%d/%d/%d %d:%d:%f %f  %f", &Y, &M, &D, &h, &m, &s, &rtasc, &decl) == EOF) {
                break;
            }
            obs[i][0] = Mjday(Y,M,D,h,m,s);
            obs[i][1]  = Rad * rtasc;
            obs[i][2]  = Rad * decl;
        }
    } else {
        printf("Es null");
    }

    fclose(fid);

    // station
    double lat = Rad*39.13607;     // [rad]
    double lon = Rad*(-121.35072); // [rad]
    double alt = 0.09981638e3;     // [m]
    double Rs[TAM];
    Position(lon, lat, alt, Rs);

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[1][0];
    double Mjd3 = obs[2][0];

    // First component
    double Mjd_UTC = Mjd1;
    double UT1_UTC, TAI_UTC, x_pole, y_pole, ddpsi, ddeps;
    IERS(eop, filas, columnas, Mjd_UTC, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    double P[ROWS][COLS], N[ROWS][COLS], polematrix[ROWS][COLS], ghamatrix[ROWS][COLS], E[ROWS][COLS];
    double aux1[ROWS][COLS], aux2[ROWS][COLS];
    PrecMatrix(MJD_J2000, Mjd_TT, P);
    NutMatrix(Mjd_TT, N);
    // E = PoleMatrix(x_pole, y_pole)*GHAMatrix(Mjd_UT1)*N*P;
    PoleMatrix(x_pole, y_pole, polematrix);
    GHAMatrix(Mjd_UT1, eop, filas, columnas, ghamatrix);
    crossMatrix(polematrix, ghamatrix, aux1);
    crossMatrix(aux1, N, aux2);
    crossMatrix(aux2, P, E);
    double Etraspuesta[ROWS][COLS];
    double rsite1[TAM];
    traspuesta(E, Etraspuesta);
    crossMatrixVector(Etraspuesta, Rs, rsite1);

    // Second component
    Mjd_UTC = Mjd2;
    IERS(eop, filas, columnas, Mjd_UTC, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    double rsite2[TAM];
    PrecMatrix(MJD_J2000, Mjd_TT, P);
    NutMatrix(Mjd_TT, N);
    // E = PoleMatrix(x_pole, y_pole)*GHAMatrix(Mjd_UT1)*N*P;
    PoleMatrix(x_pole, y_pole, polematrix);
    GHAMatrix(Mjd_UT1, eop, filas, columnas, ghamatrix);
    crossMatrix(polematrix, ghamatrix, aux1);
    crossMatrix(aux1, N, aux2);
    crossMatrix(aux2, P, E);
    traspuesta(E, Etraspuesta);
    crossMatrixVector(Etraspuesta, Rs, rsite2);

    // Third component
    Mjd_UTC = Mjd3;
    IERS(eop, filas, columnas, Mjd_UTC, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    double rsite3[TAM];
    PrecMatrix(MJD_J2000, Mjd_TT, P);
    NutMatrix(Mjd_TT, N);
    // E = PoleMatrix(x_pole, y_pole)*GHAMatrix(Mjd_UT1)*N*P;
    PoleMatrix(x_pole, y_pole, polematrix);
    GHAMatrix(Mjd_UT1, eop, filas, columnas, ghamatrix);
    crossMatrix(polematrix, ghamatrix, aux1);
    crossMatrix(aux1, N, aux2);
    crossMatrix(aux2, P, E);
    traspuesta(E, Etraspuesta);
    crossMatrixVector(Etraspuesta, Rs, rsite3);

    double *r2, *v2;
    anglesdr(obs[0][1], obs[1][1], obs[2][1], obs[0][2], obs[1][2], obs[2][2], Mjd1, Mjd2, Mjd3, rsite1, rsite2, rsite3, &r2, &v2);
    printf("------Double-R-Iteration method------\n");
    double Y_apr[6];
    Y_apr[0] = 1e-3*r2[0];
    Y_apr[1] = 1e-3*r2[1];
    Y_apr[2] = 1e-3*r2[2];
    Y_apr[3] = 1e-3*v2[0];
    Y_apr[4] = 1e-3*v2[1];
    Y_apr[5] = 1e-3*v2[2];

    for (int i = 0; i < 6; i++) {
        printf("%lf\n", Y_apr[i]);
    }
}