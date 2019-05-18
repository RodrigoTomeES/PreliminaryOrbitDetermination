//$Header$
//------------------------------------------------------------------------------
//                 PreliminaryOrbitDeterminationTest
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides a basic test for PreliminaryOrbitDetermination.c.
*
* @note
*/
//------------------------------------------------------------------------------

// Test MatLabUtilites
#include "PreliminaryOrbitDeterminationTest.h"

// cross function has a precision of 6 decimals
#define EPSILON pow(10, -6)

//Colores para los mensajes
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
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
    printf(BLUE "---- Test PreliminaryOrbitDetermination ----\n" RESET);
    // Test unit
    testUnit();

    // Test Doubler&
    testDoubler();

    // Test Angl
    testAngl();

    // Test Newtonnu
    testNewtonnu();

    // Test R_x
    testR_x();

    // Test R_y
    testR_y();

    // Test R_z
    testR_z();

    // Test Frac
    testFrac();

    // Test Timediff
    testTimediff();

    // Test Mjday
    testMjday();

    // Test Position
    testPosition();

    // Test MeanObliquity
    testMeanObliquity();

    // Test MeanObliquity
    testNutAngles();

    // Test IERS
    testIERS();

    //Test gmst
    testGmst();

    // Test Gibbs
    testGibbs();

    // Test EqnEquinox
    testEqnEquinox();

    // Test Gast
    testGast();

    // Test PoleMatrix
    testPoleMatrix();

    // Test NutMatrix
    testNutMatrix();

    // Test PrecMatrix
    testPrecMatrix();

    //Test GHAMatrix
    testGHAMatrix();

    //Test Lambert_gooding
    testLambert_gooding();

    // Test rv2coe
    testRv2coe();

    //Pasa todos los test
    printf(GREEN "---- All Pass Test From Preliminary Orbit Determination----\n" RESET);
    printf("\n");
}

void testUnit(){
    printf("---- Test UNIT ----\n");

    double v1[] = {0.0,0.0,0.0};
    double v2[]={4,-1,2};
    double v3[]={2,-2,-1};

    double v1_unitario[3];
    double v2_unitario[3];
    double v3_unitario[3];

    unit(v1, v1_unitario);
    unit(v2, v2_unitario);
    unit(v3, v3_unitario);

    //Vector v1
    printf("v1_unitario[0] = %lf, esperado 0.0\n",v1_unitario[0]);
    printf("v1_unitario[1] = %lf, esperado 0.0\n",v1_unitario[1]);
    printf("v1_unitario[2] = %lf, esperado 0.0\n",v1_unitario[2]);

    assert(fabs(v1_unitario[0] - 0.0) < EPSILON);
    assert(fabs(v1_unitario[1] - 0.0) < EPSILON);
    assert(fabs(v1_unitario[2] - 0.0) < EPSILON);

    printf("\n");

    //Vector v2
    printf("v2_unitario[0] = %lf, esperado 0.872871560943970\n",v2_unitario[0]);
    printf("v2_unitario[1] = %lf, esperado -0.218217890235992\n",v2_unitario[1]);
    printf("v2_unitario[2] = %lf, esperado 0.436435780471985\n",v2_unitario[2]);

    assert(fabs(v2_unitario[0] - 0.872871560943970) < EPSILON);
    assert(fabs(v2_unitario[1] - -0.218217890235992) < EPSILON);
    assert(fabs(v2_unitario[2] - 0.436435780471985) < EPSILON);

    printf("\n");

    //Vector v3
    printf("v3_unitario[0] = %lf, esperado 0.666666666666667\n",v3_unitario[0]);
    printf("v3_unitario[1] = %lf, esperado -0.666666666666667\n",v3_unitario[1]);
    printf("v3_unitario[2] = %lf, esperado -0.333333333333333\n",v3_unitario[2]);

    assert(fabs(v3_unitario[0] - 0.666666666666667) < EPSILON);
    assert(fabs(v3_unitario[1] - -0.666666666666667) < EPSILON);
    assert(fabs(v3_unitario[2] - -0.333333333333333) < EPSILON);

    printf(GREEN "---- Pass Test UNIT ----\n" RESET);
}

void testDoubler(){
    printf("---- Test DOUBLER ----\n");

    double r2 [TAM];
    double r3 [TAM];
    double f1;
    double f2;
    double q1;
    double magr1;
    double magr2;
    double a ;
    double deltae32;

    double cc1 = 5972180.93003294;
    double cc2 = 6395944.28126917;
    double magrsite1 = 6369760.82916145;
    double magrsite2 = 6369760.82916145;
    double magrlin = 9163883.96041412;
    double magr2in = 9163864.49341983;
    double los1[TAM] = {0.633886095165154,-0.773340379778256,0.0115358294325144};
    double los2[TAM] = {0.935539825569649,-0.0164830818224054, -0.35283642497161};
    double los3[TAM] = { 0.596384368303438,0.536072673967085,-0.597454411205648};
    double rsite1[TAM] = {4950990.3382646,256563.116260381,3999465.34658133};
    double rsite2[TAM] = { 4935037.85913036,472703.320202615, 3999475.70573182};
    double rsite3[TAM] = { 4909646.95198536,687938.936915757, 3999494.94894739};
    double t1 = -600.000004470348;
    double t3 = 600.000004470348;
    char direct = 'y';

    printf("---- Ejecutamos doubler ----\n");
    doubler(cc1,cc2,magrsite1,magrsite2,magrlin,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,
        r2,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);

    double r2Res[3]={8794373.69857176,404706.483848888,2543937.17799019};

    assert(vectoresIguales(r2,r2Res));
    double r3Res[3]={8330719.98165755,3763042.61890889,572283.772288279};
    assert(vectoresIguales(r3,r3Res));;
    assert(fabs(f1 - 0.005955295306876) < EPSILON);
    assert(fabs(f2 - -0.0256788305953251) < EPSILON);
    assert(fabs(q1 - 0.0263603467908808) < EPSILON);
    assert(fabs(magr1 - 9163883.96041412) < EPSILON);
    assert(fabs(magr2 - 9163864.49341983) < EPSILON);
    assert(fabs(a - 9138034.70407932) < EPSILON);
    assert(fabs(deltae32 - 0.432542922260443) < EPSILON);


    printf(GREEN "---- Pass Test DOUBLER ----\n" RESET);
    /*
            r2 =

                8794373.69857176
                404706.483848888
                2543937.17799019


        r3 =

                8330719.98165755
                3763042.61890889
                572283.772288279


        f1 =

                0.005955295306876


        f2 =

            -0.0256788305953251


        q1 =

                0.0263603467908808


        magr1 =

                9163883.96041412


        magr2 =

                9163864.49341983


        a =

                9138034.70407932


        deltae32 =
                        0.432542922260443
    */
}

void testAngl() {
    printf("---- Test ANGL ----\n");

    double vector1[]={2,-2,-1};
    double vector2[]={4,-1,2};

    double resultadoFuncion =  angl(vector1, vector2);
    double resultadoReal = 0.949715633622426;

    printf("  Norma función: %lf\n", resultadoFuncion);
    printf("  Norma real: %lf\n", resultadoReal);
    printf("  Diferencia: %lf\n", fabs(resultadoFuncion - resultadoReal));
    assert(fabs(resultadoFuncion - resultadoReal) < EPSILON);

    printf(GREEN "---- Pass Test ANGL ----\n" RESET);
}

void testNewtonnu() {
    printf("---- Test NEWTONNU ----\n");

    // Input
    double ecc = M_PI;
    double nu = M_PI;
    double e0, m;

    // Output
    double e0V = 9.999999000000000e+05;
    double mV = 9.999999000000000e+05;

    // Execution
    newtonnu(ecc, nu, &e0, &m);

    // Test
    printf("    --E0--\n");
    assert(fabs(e0 - e0V) < EPSILON);
    printf("  E0 función: %lf\n", e0);
    printf("  E0 real: %lf\n", e0V);
    printf("  Diferencia: %lf\n", fabs(e0 - e0V));\

    printf("\n");

    printf("    --M--\n");
    assert(fabs(m - mV) < EPSILON);
    printf("  M función: %lf\n", m);
    printf("  M real: %lf\n", mV);
    printf("  Diferencia: %lf\n", fabs(m - mV));\

    printf(GREEN "---- Pass Test NEWTONNU ----\n" RESET);
}

void testR_x() {
    printf("---- Test R_x ----\n");

    // Input
    double angle = M_PI;
    double rotmat[3][3];

    // Output
    double rotmat_result[3][3] = {{1.000000000000000,   0.000000000000000,   0.000000000000000},
                                  {0.000000000000000,  -1.000000000000000,   0.000000000000000},
                                  {0.000000000000000,  -0.000000000000000,  -1.000000000000000}};

    // Execution
    R_x(angle, rotmat);

    // Test
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("  rotmat[%d][%d] función: %lf\n", i, j, rotmat[i][j]);
            printf("  rotmat_result[%d][%d] real: %lf\n", i, j, rotmat_result[i][j]);
            printf("  Diferencia: %lf\n", fabs(rotmat[i][j] - rotmat_result[i][j]));
            assert(fabs(rotmat[i][j] - rotmat_result[i][j]) < EPSILON);
            printf("\n");
        }
    }

    printf(GREEN "---- Pass Test R_x ----\n" RESET);
}

void testR_y() {
    printf("---- Test R_y ----\n");

    // Input
    double angle = M_PI;
    double rotmat[3][3];

    // Output
    double rotmat_result[3][3] = {{-1.000000000000000,  0.000000000000000,  -0.000000000000000},
                                  { 0.000000000000000,  1.000000000000000,   0.000000000000000},
                                  { 0.000000000000000,  0.000000000000000,  -1.000000000000000}};

    // Execution
    R_y(angle, rotmat);

    // Test
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("  rotmat[%d][%d] función: %lf\n", i, j, rotmat[i][j]);
            printf("  rotmat_result[%d][%d] real: %lf\n", i, j, rotmat_result[i][j]);
            printf("  Diferencia: %lf\n", fabs(rotmat[i][j] - rotmat_result[i][j]));
            assert(fabs(rotmat[i][j] - rotmat_result[i][j]) < EPSILON);
            printf("\n");
        }
    }

    printf(GREEN "---- Pass Test R_y ----\n" RESET);
}

void testR_z() {
    printf("---- Test R_z ----\n");

    // Input
    double angle = M_PI;
    double rotmat[3][3];

    // Output
    double rotmat_result[3][3] = {{-1.000000000000000,  0.000000000000000, 0.000000000000000},
                                  {-0.000000000000000, -1.000000000000000, 0.000000000000000},
                                  { 0.000000000000000,  0.000000000000000, 1.000000000000000}};

    // Execution
    R_z(angle, rotmat);


    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("  rotmat[%d][%d] función: %lf\n", i, j, rotmat[i][j]);
            printf("  rotmat_result[%d][%d] real: %lf\n", i, j, rotmat_result[i][j]);
            printf("  Diferencia: %lf\n", fabs(rotmat[i][j] - rotmat_result[i][j]));
            assert(fabs(rotmat[i][j] - rotmat_result[i][j]) < EPSILON);
            printf("\n");
        }
    }

    printf(GREEN "---- Pass Test R_z ----\n" RESET);
}

void testFrac() {
    printf("---- Test FRAC ----\n");

    // Input
    double x = M_PI;
    double res;

    // Output
    double resultadoReal = 0.141592653589793;

    // Execution
    res=Frac(x);

    // Test
    assert(fabs(res - resultadoReal) < EPSILON);
    printf("  Frac función: %lf\n", res);
    printf("  Frac real: %lf\n", resultadoReal);
    printf("  Diferencia: %lf\n", fabs(res - resultadoReal));

    printf(GREEN "---- Pass Test FRAC ----\n" RESET);
}

void testTimediff() {
    printf("---- Test TIMEDIFF ----\n");

    // Input
    double UT1_UTC = M_PI;
    double TAI_UTC = 3*M_PI;

    // Output
    double UT1_TAI = -6.283185307179586;
    double UTC_GPS = 9.575222039230621;
    double UT1_GPS = 12.716814692820414;
    double TT_UTC  = 41.608777960769373;
    double GPS_UTC = -9.575222039230621;

    double UT1_TAIr, UTC_GPSr, UT1_GPSr, TT_UTCr, GPS_UTCr;
    // Execution
    timediff(UT1_UTC, TAI_UTC, &UT1_TAIr, &UTC_GPSr, &UT1_GPSr, &TT_UTCr, &GPS_UTCr);

    // Test
    printf("  UT1_TAI funcion: %lf\n", UT1_TAIr);
    printf("  UT1_TAI real: %lf\n", UT1_TAI);
    printf("  Diferencia: %lf\n", fabs(UT1_TAIr - UT1_TAI));
    assert(fabs(UT1_TAIr - UT1_TAI) < EPSILON);
    printf("\n");

    printf("  UTC_GPS funcion: %lf\n", UTC_GPSr);
    printf("  UTC_GPS real: %lf\n", UTC_GPS);
    printf("  Diferencia: %lf\n", fabs(UTC_GPSr - UTC_GPS));
    assert(fabs(UTC_GPSr - UTC_GPS) < EPSILON);
    printf("\n");

    printf("  UT1_GPS funcion: %lf\n", UT1_GPSr);
    printf("  UT1_GPS real: %lf\n", UT1_GPS);
    printf("  Diferencia: %lf\n", fabs(UT1_GPSr - UT1_GPS));
    assert(fabs(UT1_GPSr - UT1_GPS) < EPSILON);
    printf("\n");

    printf("  TT_UTC funcion: %lf\n", TT_UTCr);
    printf("  TT_UTC real: %lf\n", TT_UTC);
    printf("  Diferencia: %lf\n", fabs(TT_UTCr - TT_UTC));
    assert(fabs(TT_UTCr - TT_UTC) < EPSILON);
    printf("\n");

    printf("  GPS_UTC funcion: %lf\n", GPS_UTCr);
    printf("  GPS_UTC real: %lf\n", GPS_UTC);
    printf("  Diferencia: %lf\n", fabs(GPS_UTCr - GPS_UTC));
    assert(fabs(GPS_UTCr - GPS_UTC) < EPSILON);
    printf("\n");

    printf(GREEN "---- Pass Test TIMEDIFF ----\n" RESET);
}

void testMjday() {
    printf("---- Test MJDAY ----\n");

    // Input
    double year = 2015.0;
    double month = 11.0;
    double day = 3.0;
    double hour = 8.0;
    double min = 3.0;
    double sec = 4.0;

    // Output
    double Mjd_real = 5.732933546296274e+04;

    // Execution
    double Mjd = Mjday(year, month, day, hour, min, sec);

    // Test
    assert(fabs(Mjd - Mjd_real) < EPSILON);
    printf("  Mjday función: %lf\n", Mjd);
    printf("  Mjday real: %lf\n", Mjd_real);
    printf("  Diferencia: %lf\n", fabs(Mjd - Mjd_real));

    printf(GREEN "---- Pass Test MJDAY ----\n" RESET);
}

void testPosition() {
    printf("---- Test POSITION ----\n");

    // Input
    double lon = M_PI;
    double lat = M_PI;
    double h = 40;
    double r[3];

    // Output
    double r_result[] = {6378177.000000000, -7.811014047445268e-10, 7.758724479233378e-10};

    // Execution
    Position(lon, lat, h, r);

    // Test
    for (int i = 0; i < 3; i++) {
        printf("  r[%d] función: %lf\n", i, r[i]);
        printf("  r_result[%d] real: %lf\n", i, r_result[i]);
        printf("  Diferencia: %lf\n", fabs(r[i] - r_result[i]));
        assert(fabs(r[i] - r_result[i]) < EPSILON);
        printf("\n");
    }

    printf(GREEN "---- Pass Test POSITION ----\n" RESET);
}

void testMeanObliquity() {
    printf("---- Test MEAN OBLIQUITY ----\n");

    // Input
    double Mjd_TT = 5.732933546296274e+04;
    double MOblq;

    // Output
    double MOblq_real = 0.409056857329232;

    // Execution
    MOblq = MeanObliquity(Mjd_TT);

    // Test
    printf("  MeanObliquity función: %lf\n", MOblq);
    printf("  MeanObliquity real: %lf\n", MOblq_real);
    printf("  Diferencia: %lf\n", fabs(MOblq - MOblq_real));
    assert(fabs(MOblq - MOblq_real) < EPSILON);

    printf(GREEN "---- Pass Test MEAN OBLIQUITY ----\n" RESET);
}

void testNutAngles() {
    printf("---- Test NUT ANGLES ----\n");

    // Input
    double Mjd_TT = 5.732933546296274e+04;
    double dpsi,deps;

    // Output
    double NutAngles_const_real[] = {-7.840339161497435e-06, -4.469852616699864e-05};

    // Execution
    NutAngles(Mjd_TT, &dpsi, &deps);

    // Test
    printf("  dpsi función: %lf\n", dpsi);
    printf("  dpsi real: %lf\n", NutAngles_const_real[0]);
    printf("  Diferencia: %lf\n", fabs(dpsi - NutAngles_const_real[0]));
    assert(fabs(dpsi - NutAngles_const_real[0]) < EPSILON);
    printf("\n");

    printf("  deps función: %lf\n", deps);
    printf("  deps real: %lf\n", NutAngles_const_real[1]);
    printf("  Diferencia: %lf\n", fabs(deps - NutAngles_const_real[1]));
    assert(fabs(deps - NutAngles_const_real[1]) < EPSILON);

    printf(GREEN "---- Pass Test NUT ANGLES ----\n" RESET);
}

void testIERS(){
    printf("---- Test IERS ----\n");

    FILE* fid = fopen("eop19620101.txt","rt");

    int filas = 20026;
    int columnas = 13;
    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;

    if (fid == NULL){
        exit(EXIT_FAILURE);
    }

    double **eop;
    eop = (double **) malloc (filas*sizeof(double *));

    if (eop != NULL) {
        for (int i = 0; i < filas; i++) {
            eop[i] = (double *) malloc (columnas * sizeof(double));
            if (fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) == EOF) {
                break;
            }
            eop[i][0]  =  v1;
            eop[i][1]  =  v2;
            eop[i][2]  =  v3;
            eop[i][3]  =  v4;
            eop[i][4]  =  v5;
            eop[i][5]  =  v6;
            eop[i][6]  =  v7;
            eop[i][7]  =  v8;
            eop[i][8]  =  v9;
            eop[i][9]  = v10;
            eop[i][10] = v11;
            eop[i][11] = v12;
            eop[i][12] = v13;

            //printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }
    } else {
        printf("Es null");
    }

    fclose(fid);

    double ** traspuesta;
    traspuesta=(double **) malloc (columnas*sizeof(double *));
    double traspuestaFilas = columnas;
    double traspuestaColumnas = filas;

    if (traspuesta != NULL) {
        for (int i = 0; i < columnas; i++) {
            traspuesta[i] = (double *) malloc (filas * sizeof(double));
            for(int j=0;j<filas;j++){
                traspuesta[i][j]=eop[j][i];
            }
        }
    }

    // for(int j=0;j<filas;j++){
    //     for (int i = 0; i < columnas; i++) {
    //         printf("%lf ", traspuesta[i][j]);
    //     }
    //     printf("\n");
    // }

    // Input
    // double **eop;
    double Mjd_UTC = 54977.6669036457;
    char interp = 'l';
    double UT1_UTC;
    double TAI_UTC;
    double x_pole;
    double y_pole;
    double ddpsi;
    double ddeps;

    // Output
    double UT1_UTC_real = 0.258022690875596;
    double TAI_UTC_real = 34;
    double x_pole_real = 7.5789200806793e-08;
    double y_pole_real = 2.56777193042581e-06;
    double ddpsi_real = -2.91335538497448e-07;
    double ddeps_real = -4.54223178651006e-08;

    // Execution
    IERS(traspuesta, traspuestaFilas, traspuestaColumnas,Mjd_UTC, interp, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);

    // Test
    printf("  UT1_UTC función: %lf\n", UT1_UTC);
    printf("  UT1_UTC real: %lf\n", UT1_UTC_real);
    printf("  Diferencia: %lf\n", fabs(UT1_UTC - UT1_UTC_real));
    assert(fabs(UT1_UTC - UT1_UTC_real) < EPSILON);
    printf("\n");

    printf("  TAI_UTC función: %lf\n", TAI_UTC);
    printf("  TAI_UTC real: %lf\n", TAI_UTC_real);
    printf("  Diferencia: %lf\n", fabs(TAI_UTC - TAI_UTC_real));
    assert(fabs(TAI_UTC - TAI_UTC_real) < EPSILON);
    printf("\n");

    printf("  x_pole función: %lf\n", x_pole);
    printf("  x_pole real: %lf\n", x_pole_real);
    printf("  Diferencia: %lf\n", fabs(x_pole - x_pole_real));
    assert(fabs(x_pole - x_pole_real) < EPSILON);
    printf("\n");

    printf("  y_pole función: %lf\n", y_pole);
    printf("  y_pole real: %lf\n", y_pole_real);
    printf("  Diferencia: %lf\n", fabs(y_pole - y_pole_real));
    assert(fabs(y_pole - y_pole_real) < EPSILON);
    printf("\n");

    printf("  ddpsi función: %lf\n", ddpsi);
    printf("  ddpsi real: %lf\n", ddpsi_real);
    printf("  Diferencia: %lf\n", fabs(ddpsi - ddpsi_real));
    assert(fabs(ddpsi - ddpsi_real) < EPSILON);
    printf("\n");

    printf("  ddeps función: %lf\n", ddeps);
    printf("  ddeps real: %lf\n", ddeps_real);
    printf("  Diferencia: %lf\n", fabs(ddeps - ddeps_real));
    assert(fabs(ddeps - ddeps_real) < EPSILON);

    printf(GREEN "---- Pass Test IERS ----\n" RESET);
}

void testGmst(){

    printf("---- Test gmst ----\n");

    // Input
    double Mjd_UT1 = 54977.680795521;
    double gmstV;

    // Output
    double gmst_real = 2.25937441720215;

    // Execution
    gmstV = gmst(Mjd_UT1);

    // Test
    printf("  gmst función: %lf\n", gmstV);
    printf("  gmst real: %lf\n", gmst_real);
    printf("  Diferencia: %lf\n", fabs(gmstV - gmst_real));
    assert(fabs(gmstV - gmst_real) < EPSILON);

    printf(GREEN "---- Pass Test gmst ----\n" RESET);
}

void testGibbs() {
    printf("---- Test GIBBS ----\n");

    // Input
    double r1[] = {20387627.0717529,1865163.69633398,-109943.688555879};
    double r2[] = {20435422.3521544,1070699.44671825,1012905.49143365};
    double r3[] = {20398157.0666256,271778.615869788,2131538.39542076};
    double v2[3];
    double theta;
    double theta1;
    double copa;
    char error[20];

    // Output
    double v2_result[] = {17.4448460090308,-2659.68695020331,3741.47770465728};
    double theta_result = 0.0672088229314286;
    double theta1_result = 0.0670842334897834;
    double copa_result = -7.87130777224476e-16;
    char error_result[] = "ok";

    // Execution
    gibbs(r1, r2, r3, v2, &theta, &theta1, &copa, error);

    // Test
    printf("  v2\n");
    muestraVector(v2);
    printf("  v2_result\n");
    muestraVector(v2_result);
    printf("  v2 iguales: %d\n", vectoresIguales(v2,v2_result));
    assert(vectoresIguales(v2,v2_result));
    printf("\n");

    printf("  theta función: %lf\n", theta);
    printf("  theta real: %lf\n", theta_result);
    printf("  Diferencia: %lf\n", fabs(theta - theta_result));
    assert(fabs(theta - theta_result) < EPSILON);
    printf("\n");

    printf("  theta1 función: %lf\n", theta1);
    printf("  theta1 real: %lf\n", theta1_result);
    printf("  Diferencia: %lf\n", fabs(theta1 - theta1_result));
    assert(fabs(theta1 - theta1_result) < EPSILON);
    printf("\n");

    printf("  copa función: %lf\n", copa);
    printf("  copa real: %lf\n", copa_result);
    printf("  Diferencia: %lf\n", fabs(copa - copa_result));
    assert(fabs(copa - copa_result) < EPSILON);
    printf("\n");

    printf("  error igual: %d\n", vectoresIguales(v2,v2_result));
    assert(vectoresIguales(v2,v2_result));

    printf(GREEN "---- Pass Test GIBBS ----\n" RESET);
}

void testEqnEquinox(){
    printf("---- Test EqnEquinox ----\n");

    // Input
    double Mjd_TT = 54977.6815585532;
    double EqnEquinoxV;

    // Output
    double EqnEquinox_real = 5.95287721905946e-05;

    // Execution
    EqnEquinoxV = EqnEquinox(Mjd_TT);

    // Test
    printf("  EqnEquinox función: %lf\n", EqnEquinoxV);
    printf("  EqnEquinox real: %lf\n", EqnEquinox_real);
    printf("  Diferencia: %lf\n", fabs(EqnEquinoxV - EqnEquinox_real));
    assert(fabs(EqnEquinoxV - EqnEquinox_real) < EPSILON);

    printf(GREEN "---- Pass Test EqnEquinox ----\n" RESET);
}

void testGast() {
    printf("---- Test GAST ----\n");

    FILE* fid = fopen("eop19620101.txt","rt");

    int filas = 20026;
    int columnas = 13;
    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;

    if (fid == NULL){
        exit(EXIT_FAILURE);
    }

    double **eop;
    eop = (double **) malloc (filas*sizeof(double *));

    if (eop != NULL) {
        for (int i = 0; i < filas; i++) {
            eop[i] = (double *) malloc (columnas * sizeof(double));
            if (fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) == EOF) {
                break;
            }
            eop[i][0]  =  v1;
            eop[i][1]  =  v2;
            eop[i][2]  =  v3;
            eop[i][3]  =  v4;
            eop[i][4]  =  v5;
            eop[i][5]  =  v6;
            eop[i][6]  =  v7;
            eop[i][7]  =  v8;
            eop[i][8]  =  v9;
            eop[i][9]  = v10;
            eop[i][10] = v11;
            eop[i][11] = v12;
            eop[i][12] = v13;

            //printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }
    } else {
        printf("Es null");
    }

    fclose(fid);

    double ** traspuesta;
    traspuesta=(double **) malloc (columnas*sizeof(double *));
    double traspuestaFilas = columnas;
    double traspuestaColumnas = filas;

    if (traspuesta != NULL) {
        for (int i = 0; i < columnas; i++) {
            traspuesta[i] = (double *) malloc (filas * sizeof(double));
            for(int j=0;j<filas;j++){
                traspuesta[i][j]=eop[j][i];
            }
        }
    }

    // for(int j=0;j<filas;j++){
    //     for (int i = 0; i < columnas; i++) {
    //         printf("%lf ", traspuesta[i][j]);
    //     }
    //     printf("\n");
    // }

    // Input
    // double **eop;
    double Mjd_UT1 = 54977.680795521;
    double gstime;

    // Output
    double gstime_real = 2.25943394597434;

    // Execution
    gstime = gast(Mjd_UT1, traspuesta, traspuestaFilas, traspuestaColumnas);

    // Test
    printf("  gstime función: %lf\n", gstime);
    printf("  gstime real: %lf\n", gstime_real);
    printf("  Diferencia: %lf\n", fabs(gstime - gstime_real));
    assert(fabs(gstime - gstime_real) < EPSILON);

    printf(GREEN "---- Pass Test GAST ----\n" RESET);
}

void testPoleMatrix(){
    printf("---- Test PoleMatrix ----\n");

    // Input

    double PoleMat[3][3];

    // Output
    double PoleMat_result[3][3] = {{0.999999999999997,  1.94609382460874e-13, 7.57892008065431e-08},
                                  {0, 0.999999999996703, -2.56777193042298e-06},
                                  { -7.57892008067929e-08,  2.56777193042298e-06, 0.9999999999967}};

    // Execution
    PoleMatrix(7.5789200806793e-08,2.56777193042581e-06,PoleMat);


    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("  PoleMat[%d][%d] función: %lf\n", i, j, PoleMat[i][j]);
            printf("  PoleMat_result[%d][%d] real: %lf\n", i, j, PoleMat_result[i][j]);
            printf("  Diferencia: %lf\n", fabs(PoleMat[i][j] - PoleMat_result[i][j]));
            assert(fabs(PoleMat[i][j] - PoleMat_result[i][j]) < EPSILON);
            printf("\n");
        }
    }

    printf(GREEN "---- Pass Test PoleMatrix ----\n" RESET);
}

void testNutMatrix(){
    printf("---- Test NUT MATRIX ----\n");

    // Input
    double Mjd_TT = 54977.6676696643;
    double NutMat[3][3];

    // Output
    double NutMat_real[3][3] = {{   0.999999997895984,     -5.95170051004573e-05,     -2.58022713429399e-05},
                                {5.95164295625406e-05,         0.999999997980121,     -2.23059014984558e-05},
                                {2.58035988712757e-05,      2.23043657924436e-05,         0.999999999418345}};

    // Execution
    NutMatrix(Mjd_TT, NutMat);

    // Test
    muestraMatriz(NutMat);
    printf("  NutMat iguales: %d\n", matricesIguales(NutMat, NutMat_real));
    assert(matricesIguales(NutMat, NutMat_real));

    printf(GREEN "---- Pass Test NUT MATRIX ----\n" RESET);
}

void testPrecMatrix(){
    printf("---- Test PrecMatrix ----\n");

    // Input

    double PrecMat[3][3];

    // Output
    double PrecMat_result[3][3] = {{0.99999737378108,       -0.0021019566973475,     -0.000913350417381566},
                                  {0.00210195669733333,          0.99999779088612,     -9.59928282397388e-07},
                                  { 0.000913350417414164,     -9.59897265411807e-07,          0.99999958289496}};

    // Execution
    PrecMatrix(51544.5,54977.6815585532,PrecMat);


    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("  PrecMat[%d][%d] función: %lf\n", i, j, PrecMat[i][j]);
            printf("  PrecMat_result[%d][%d] real: %lf\n", i, j, PrecMat_result[i][j]);
            printf("  Diferencia: %lf\n", fabs(PrecMat[i][j] - PrecMat_result[i][j]));
            assert(fabs(PrecMat[i][j] - PrecMat_result[i][j]) < EPSILON);
            printf("\n");
        }
    }

    printf(GREEN "---- Pass Test PrecMatrix ----\n" RESET);
}

void testGHAMatrix(){
    printf("---- Test GHA_MATRIX ----\n");

    FILE* fid = fopen("eop19620101.txt","rt");

    int filas = 20026;
    int columnas = 13;
    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;

    if (fid == NULL){
        exit(EXIT_FAILURE);
    }

    double **eop;
    eop = (double **) malloc (filas*sizeof(double *));

    if (eop != NULL) {
        for (int i = 0; i < filas; i++) {
            eop[i] = (double *) malloc (columnas * sizeof(double));
            if (fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) == EOF) {
                break;
            }
            eop[i][0]  =  v1;
            eop[i][1]  =  v2;
            eop[i][2]  =  v3;
            eop[i][3]  =  v4;
            eop[i][4]  =  v5;
            eop[i][5]  =  v6;
            eop[i][6]  =  v7;
            eop[i][7]  =  v8;
            eop[i][8]  =  v9;
            eop[i][9]  = v10;
            eop[i][10] = v11;
            eop[i][11] = v12;
            eop[i][12] = v13;

            //printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }
    } else {
        printf("Es null");
    }

    fclose(fid);

    double ** traspuesta;
    traspuesta=(double **) malloc (columnas*sizeof(double *));
    double traspuestaFilas = columnas;
    double traspuestaColumnas = filas;

    if (traspuesta != NULL) {
        for (int i = 0; i < columnas; i++) {
            traspuesta[i] = (double *) malloc (filas * sizeof(double));
            for(int j=0;j<filas;j++){
                traspuesta[i][j]=eop[j][i];
            }
        }
    }

    // for(int j=0;j<filas;j++){
    //     for (int i = 0; i < columnas; i++) {
    //         printf("%f ", traspuesta[i][j]);
    //     }
    //     printf("\n");
    // }

    // Input
    double Mjd_UT1 = 54977.680795521;
    double GHAmat[3][3];

    // Output
    double GHAmat_real[3][3] = {{-0.635485861035896,         0.772112505029847,                         0},
                                {-0.772112505029847,        -0.635485861035896,                         0},
                                {                 0,                         0,                         1}};

    // Execution
    GHAMatrix(Mjd_UT1, traspuesta, traspuestaFilas, traspuestaColumnas, GHAmat);

    // Test
    printf("  GHAmat\n");
    muestraMatriz(GHAmat);
    printf("  GHAmat_real\n");
    muestraMatriz(GHAmat_real);
    printf("  GHAmat iguales: %d\n", matricesIguales(GHAmat,GHAmat_real));
    assert(matricesIguales(GHAmat,GHAmat_real));

	printf(GREEN "---- Pass Test GHA_MATRIX ----\n" RESET);
}

void testLambert_gooding(){
    testTLamb();
    testd8rt();
    testXLamb();
    testVLamb();
}

void testTLamb(){
    printf("---- Test tlamb ----\n");

    // Input
    double t;
    double dt;
    double d2t;
    double d3t;

    // Output
    double t_result = 8.22665277455477;
    double dt_result = -3.54465085248612e-11;
    double d2t_result =  28.5802140304226;
    double d3t_result = 26.4789138610536;

    // Execution
    tlamb(1,0.804611564466232,0.352600230327203,0.142219617787892,3, &t, &dt, &d2t, &d3t);

    // Test
    printf("  t función: %lf\n", t);
    printf("  t real: %lf\n", t_result);
    printf("  Diferencia: %lf\n", fabs(t - t_result));
    assert(fabs(t - t_result) < EPSILON);
    printf("\n");

    printf("  dt función: %lf\n", dt);
    printf("  dt real: %lf\n", dt_result);
    printf("  Diferencia: %lf\n", fabs(dt - dt_result));
    assert(fabs(dt - dt_result) < EPSILON);
    printf("\n");

    printf("  d2t función: %lf\n", d2t);
    printf("  d2t real: %lf\n", d2t_result);
    printf("  Diferencia: %lf\n", fabs(d2t - d2t_result));
    assert(fabs(d2t - d2t_result) < EPSILON);
    printf("\n");

    printf("  d3t función: %lf\n", d3t);
    printf("  d3t real: %lf\n", d3t_result);
    printf("  Diferencia: %lf\n", fabs(d3t - d3t_result));
    assert(fabs(d3t - d3t_result) < EPSILON);

    printf(GREEN "---- Pass Test tlamb ----\n" RESET);
}

void testd8rt(){
    printf("---- Test d8rt ----\n");

    // Input
    double x;

    // Output
    double x_result = 0.780220027563195;

    // Execution
    x = d8rt(0.137320935252476);

    // Test
    printf("  t función: %lf\n", x);
    printf("  t real: %lf\n", x_result);
    printf("  Diferencia: %lf\n", fabs(x - x_result));
    assert(fabs(x - x_result) < EPSILON);

    printf(GREEN "---- Pass Test d8rt ----\n" RESET);
}

void testXLamb(){
    printf("---- Test xlamb ----\n");

    // Input
    double n;
    double x1;
    double x2;

    // Output
    double n_result = 0;
    double x1_result = 0;
    double x2_result =  0;

    // Execution
    xlamb(1,0.804611564466232,0.352600230327203,0.913438160983316, &n, &x1, &x2);

    // Test
    printf("  n función: %lf\n", n);
    printf("  n real: %lf\n", n_result);
    printf("  Diferencia: %lf\n", fabs(n - n_result));
    assert(fabs(n - n_result) < EPSILON);
    printf("\n");

    printf("  x1 función: %lf\n", x1);
    printf("  x1 real: %lf\n", x1_result);
    printf("  Diferencia: %lf\n", fabs(x1 - x1_result));
    assert(fabs(x1 - x1_result) < EPSILON);
    printf("\n");

    printf("  x2 función: %lf\n", x2);
    printf("  x2 real: %lf\n", x2_result);
    printf("  Diferencia: %lf\n", fabs(x2 - x2_result));
    assert(fabs(x2 - x2_result) < EPSILON);
    printf("\n");

    printf(GREEN "---- Pass Test xlamb ----\n" RESET);
}

void testVLamb(){
    printf("---- Test vlamb ----\n");

    // Input

    double n;
    double * vri;
    double * vti;
    double * vrf;
    double * vtf;
    


    // Output
    double n_result = 1;
    double vri_result[] = { -4.24634165279188,0};
    double vti_result[] = { 6585.8979308205,0};
    double vrf_result[] = { -11.6872756516642,0};
    double vtf_result[] = { 6589.38917620768,0};
    

    double r1[]={9163781.53546157,0,0};
    double r2[]={9158926.30394324,0,0};

    vlamb(398600441800000,r1,r2,0.431406120912815,600.000004470348, &n, &vri,&vti,&vrf,&vtf);

    // Test


    printf("  n función: %f\n", n);
    printf("  n real: %f\n", n_result);
    printf("  Diferencia: %f\n", fabs(n - n_result));
    assert(fabs(n - n_result) < EPSILON);
    printf("\n");

    for(int aux=0;aux<n;aux++){
        printf("  vri función: %f\n", vri[aux]);
        printf("  vri real: %f\n", vri_result[aux]);
        printf("  Diferencia: %f\n", fabs(vri[aux] - vri_result[aux]));
        assert(fabs(vri[aux] - vri_result[aux]) < EPSILON);
        printf("\n");
    }

    for(int aux=0;aux<n;aux++){
        printf("  vti función: %f\n", vti[aux]);
        printf("  vti real: %f\n", vti_result[aux]);
        printf("  Diferencia: %f\n", fabs(vti[aux] - vti_result[aux]));
        assert(fabs(vti[aux] - vti_result[aux]) < EPSILON);
        printf("\n");
    }

    for(int aux=0;aux<n;aux++){
        printf("  vrf función: %f\n", vrf[aux]);
        printf("  vrf real: %f\n", vrf_result[aux]);
        printf("  Diferencia: %f\n", fabs(vrf[aux] - vrf_result[aux]));
        assert(fabs(vrf[aux] - vrf_result[aux]) < EPSILON);
        printf("\n");
    }

    for(int aux=0;aux<n;aux++){
        printf("  vtf función: %f\n", vtf[aux]);
        printf("  vtf real: %f\n", vtf_result[aux]);
        printf("  Diferencia: %f\n", fabs(vtf[aux] - vtf_result[aux]));
        assert(fabs(vtf[aux] - vtf_result[aux]) < EPSILON);
        printf("\n");
    }

    printf(GREEN "---- Pass Test vlamb ----\n" RESET);

}
void testRv2coe() {
    printf("---- Test RV2COE ----\n");

    // Input
    double r[] = {20418280.3742389, 1067836.39923722, 1015404.95114477};
    double v[] = {16.8797950290867, -2654.08002932654, 3734.1200461541};

    // double r[] = {3,4,5};
    // double v[] = {10,6,8};
    double p;
    double a;
    double ecc;
    double incl;
    double omega;
    double argp;
    double nu;
    double m;
    double arglat;
    double truelon;
    double lonper;

    // Output
    double p_real = 22062031.7359055;
    double a_real = 22201041.6868316;
    double ecc_real = 0.0791291077785513;
    double incl_real = 2.18716682303266;
    double omega_real = 0.0874406813869261;
    double argp_real = 6.15374268499539;
    double nu_real = 0.190267237474694;
    double m_real = 0.16199787056862;
    double arglat_real = 999999.1;
    double truelon_real = 999999.1;
    double lonper_real = 999999.1;

    // double p_real = 2.92021753599572e-12;
    // double a_real = 3.53553390593901;
    // double ecc_real = 0.999999999999587;
    // double incl_real = 2.27159902873729;
    // double omega_real = 3.06482076232001;
    // double argp_real = 5.10178346758267;
    // double nu_real = 3.14159265358979;
    // double m_real = 5.82562114015612;
    // double arglat_real = 999999.1;
    // double truelon_real = 999999.1;
    // double lonper_real = 999999.1;

    // Execution
    rv2coe(r, v, &p, &a, &ecc, &incl, &omega, &argp, &nu, &m, &arglat, &truelon, &lonper);

    // Test
    printf("  p función: %lf\n", p);
    printf("  p real: %lf\n", p_real);
    printf("  Diferencia: %lf\n", fabs(p - p_real));
    assert(fabs(p - p_real) < EPSILON);
    printf("\n");

    printf("  a función: %.15f\n", a);
    printf("  a real: %.15f\n", a_real);
    printf("  Diferencia: %lf\n", fabs(a - a_real));
    assert(fabs(a - a_real) < EPSILON);
    printf("\n");

    printf("  ecc función: %lf\n", ecc);
    printf("  ecc real: %lf\n", ecc_real);
    printf("  Diferencia: %lf\n", fabs(ecc - ecc_real));
    assert(fabs(ecc - ecc_real) < EPSILON);
    printf("\n");

    printf("  incl función: %lf\n", incl);
    printf("  incl real: %lf\n", incl_real);
    printf("  Diferencia: %lf\n", fabs(incl - incl_real));
    assert(fabs(incl - incl_real) < EPSILON);
    printf("\n");

    printf("  omega función: %lf\n", omega);
    printf("  omega real: %lf\n", omega_real);
    printf("  Diferencia: %lf\n", fabs(omega - omega_real));
    assert(fabs(omega - omega_real) < EPSILON);
    printf("\n");

    printf("  argp función: %lf\n", argp);
    printf("  argp real: %lf\n", argp_real);
    printf("  Diferencia: %lf\n", fabs(argp - argp_real));
    assert(fabs(argp - argp_real) < EPSILON);
    printf("\n");

    printf("  nu función: %lf\n", nu);
    printf("  nu real: %lf\n", nu_real);
    printf("  Diferencia: %lf\n", fabs(nu - nu_real));
    assert(fabs(nu - nu_real) < EPSILON);
    printf("\n");

    printf("  m función: %lf\n", m);
    printf("  m real: %lf\n", m_real);
    printf("  Diferencia: %lf\n", fabs(m - m_real));
    assert(fabs(m - m_real) < EPSILON);
    printf("\n");

    printf("  arglat función: %lf\n", arglat);
    printf("  arglat real: %lf\n", arglat_real);
    printf("  Diferencia: %lf\n", fabs(arglat - arglat_real));
    assert(fabs(arglat - arglat_real) < EPSILON);
    printf("\n");

    printf("  truelon función: %lf\n", truelon);
    printf("  truelon real: %lf\n", truelon_real);
    printf("  Diferencia: %lf\n", fabs(truelon - truelon_real));
    assert(fabs(truelon - truelon_real) < EPSILON);
    printf("\n");

    printf("  lonper función: %lf\n", lonper);
    printf("  lonper real: %lf\n", lonper_real);
    printf("  Diferencia: %lf\n", fabs(lonper - lonper_real));
    assert(fabs(lonper - lonper_real) < EPSILON);

    printf(GREEN "---- Pass Test RV2COE ----\n" RESET);
}
