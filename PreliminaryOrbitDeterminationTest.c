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

#define EPSILON pow(10, -7)

// Epsilon es de 10⁻7 ya que es la precisión que deja ver matlab para el test

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

    printf("  Norma función: %f\n", resultadoFuncion);
    printf("  Norma real: %f\n", resultadoReal);
    printf("  Diferencia: %f\n", fabs(resultadoFuncion - resultadoReal));
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
    printf("  E0 función: %f\n", e0);
    printf("  E0 real: %f\n", e0V);
    printf("  Diferencia: %f\n", fabs(e0 - e0V));\

    printf("\n");

    printf("    --M--\n");
    assert(fabs(m - mV) < EPSILON);
    printf("  M función: %f\n", m);
    printf("  M real: %f\n", mV);
    printf("  Diferencia: %f\n", fabs(m - mV));\

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
            printf("  rotmat[%d][%d] función: %f\n", i, j, rotmat[i][j]);
            printf("  rotmat_result[%d][%d] real: %f\n", i, j, rotmat_result[i][j]);
            printf("  Diferencia: %f\n", fabs(rotmat[i][j] - rotmat_result[i][j]));
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
            printf("  rotmat[%d][%d] función: %f\n", i, j, rotmat[i][j]);
            printf("  rotmat_result[%d][%d] real: %f\n", i, j, rotmat_result[i][j]);
            printf("  Diferencia: %f\n", fabs(rotmat[i][j] - rotmat_result[i][j]));
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
            printf("  rotmat[%d][%d] función: %f\n", i, j, rotmat[i][j]);
            printf("  rotmat_result[%d][%d] real: %f\n", i, j, rotmat_result[i][j]);
            printf("  Diferencia: %f\n", fabs(rotmat[i][j] - rotmat_result[i][j]));
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
    printf("  Frac función: %f\n", res);
    printf("  Frac real: %f\n", resultadoReal);
    printf("  Diferencia: %f\n", fabs(res - resultadoReal));

    printf(GREEN "---- Pass Test FRAC ----\n" RESET);
}

void testTimediff() {
    printf("---- Test TIMEDIFF ----\n");

    // Input
    double UT1_UTC = M_PI;
    double TAI_UTC = 3*M_PI;
    double timediff_const[5];

    // Output
    double UT1_TAI = -6.283185307179586;
    double UTC_GPS = 9.575222039230621;
    double UT1_GPS = 12.716814692820414;
    double TT_UTC  = 41.608777960769373;
    double GPS_UTC = -9.575222039230621;

    // Execution
    timediff(UT1_UTC, TAI_UTC, timediff_const);

    // Test
    printf("  UT1_TAI funcion: %f\n", timediff_const[0]);
    printf("  UT1_TAI real: %f\n", UT1_TAI);
    printf("  Diferencia: %f\n", fabs(timediff_const[0] - UT1_TAI));
    assert(fabs(timediff_const[0] - UT1_TAI) < EPSILON);
    printf("\n");

    printf("  UTC_GPS funcion: %f\n", timediff_const[1]);
    printf("  UTC_GPS real: %f\n", UTC_GPS);
    printf("  Diferencia: %f\n", fabs(timediff_const[1] - UTC_GPS));
    assert(fabs(timediff_const[1] - UTC_GPS) < EPSILON);
    printf("\n");

    printf("  UT1_GPS funcion: %f\n", timediff_const[2]);
    printf("  UT1_GPS real: %f\n", UT1_GPS);
    printf("  Diferencia: %f\n", fabs(timediff_const[2] - UT1_GPS));
    assert(fabs(timediff_const[2] - UT1_GPS) < EPSILON);
    printf("\n");

    printf("  TT_UTC funcion: %f\n", timediff_const[3]);
    printf("  TT_UTC real: %f\n", TT_UTC);
    printf("  Diferencia: %f\n", fabs(timediff_const[3] - TT_UTC));
    assert(fabs(timediff_const[3] - TT_UTC) < EPSILON);
    printf("\n");

    printf("  GPS_UTC funcion: %f\n", timediff_const[4]);
    printf("  GPS_UTC real: %f\n", GPS_UTC);
    printf("  Diferencia: %f\n", fabs(timediff_const[4] - GPS_UTC));
    assert(fabs(timediff_const[4] - GPS_UTC) < EPSILON);
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
    printf("  Mjday función: %f\n", Mjd);
    printf("  Mjday real: %f\n", Mjd_real);
    printf("  Diferencia: %f\n", fabs(Mjd - Mjd_real));

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
        printf("  r[%d] función: %f\n", i, r[i]);
        printf("  r_result[%d] real: %f\n", i, r_result[i]);
        printf("  Diferencia: %f\n", fabs(r[i] - r_result[i]));
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
    printf("  MeanObliquity función: %f\n", MOblq);
    printf("  MeanObliquity real: %f\n", MOblq_real);
    printf("  Diferencia: %f\n", fabs(MOblq - MOblq_real));
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
    printf("  dpsi función: %f\n", dpsi);
    printf("  dpsi real: %f\n", NutAngles_const_real[0]);
    printf("  Diferencia: %f\n", fabs(dpsi - NutAngles_const_real[0]));
    assert(fabs(dpsi - NutAngles_const_real[0]) < EPSILON);
    printf("\n");

    printf("  deps función: %f\n", deps);
    printf("  deps real: %f\n", NutAngles_const_real[1]);
    printf("  Diferencia: %f\n", fabs(deps - NutAngles_const_real[1]));
    assert(fabs(deps - NutAngles_const_real[1]) < EPSILON);

    printf(GREEN "---- Pass Test NUT ANGLES ----\n" RESET);
}

void testIERS(){

}

void testGmst(){

    printf("---- Test gmst ----\n");

    // Input
    double Mjd_UT1 = 54977.6738510765;
    double gmstV;

    // Output
    double gmst_real = 2.21562172211082;

    // Execution
    gmstV = gmst(Mjd_UT1);

    // Test
    printf("  gmst función: %f\n", gmstV);
    printf("  gmst real: %f\n", gmst_real);
    printf("  Diferencia: %f\n", fabs(gmstV - gmst_real));
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

    printf("  theta función: %f\n", theta);
    printf("  theta real: %f\n", theta_result);
    printf("  Diferencia: %f\n", fabs(theta - theta_result));
    assert(fabs(theta - theta_result) < EPSILON);
    printf("\n");

    printf("  theta1 función: %f\n", theta1);
    printf("  theta1 real: %f\n", theta1_result);
    printf("  Diferencia: %f\n", fabs(theta1 - theta1_result));
    assert(fabs(theta1 - theta1_result) < EPSILON);
    printf("\n");

    printf("  copa función: %f\n", copa);
    printf("  copa real: %f\n", copa_result);
    printf("  Diferencia: %f\n", fabs(copa - copa_result));
    assert(fabs(copa - copa_result) < EPSILON);
    printf("\n");

    printf("  error igual: %d\n", vectoresIguales(v2,v2_result));
    assert(vectoresIguales(v2,v2_result));
    printf("\n");

    printf(GREEN "---- Pass Test GIBBS ----\n" RESET);
}
