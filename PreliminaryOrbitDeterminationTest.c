//$Header$ 
//------------------------------------------------------------------------------ 
//                           PreliminaryOrbitDeterminationTest
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination. 
// 
// Legal: MIT  License
// 
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27 
// 
/**  
* Provides a basic test for all function of Preliminary Orbit Determination Proyect.  
*  
* @note     
*/ 
//------------------------------------------------------------------------------

// Test MatLabUtilites
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "doubler.h"

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

double v1[] = {0.0,0.0,0.0};
double v2[]={4,-1,2};
double v3[]={2,-2,-1};

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

    // Test Doubler
    testDoubler();

}

void testUnit(){
    printf("---- Test UNIT ----\n");

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


    //Vector v2
    printf("v2_unitario[0] = %lf, esperado 0.872871560943970\n",v2_unitario[0]);
    printf("v2_unitario[1] = %lf, esperado -0.218217890235992\n",v2_unitario[1]);
    printf("v2_unitario[2] = %lf, esperado 0.436435780471985\n",v2_unitario[2]);

    assert(fabs(v2_unitario[0] - 0.872871560943970) < EPSILON);
    assert(fabs(v2_unitario[1] - -0.218217890235992) < EPSILON);
    assert(fabs(v2_unitario[2] - 0.436435780471985) < EPSILON);

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
    double r2 [TAM];
    double r3 [TAM];
    double * f1;
    double *f2;
    double *q1;
    double * magr1;
    double * magr2;
    double *a ;
    double * deltae32;

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
    char direct = "y";

    doubler(cc1,cc2,magrsite1,magrsite2,magrlin,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct,
        r2,r3,f1,f2,q1,magr1,magr2,a,deltae32);

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