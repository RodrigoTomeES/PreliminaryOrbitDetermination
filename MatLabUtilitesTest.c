//$Header$ 
//------------------------------------------------------------------------------ 
//                           MatLabUtilitesTest
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination. 
// 
// Legal: MIT  License
// 
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27 
// 
/**  
* Provides a basic test for MatLabUtilites.c.  
*  
* @note     
*/ 
//------------------------------------------------------------------------------

// Test MatLabUtilites
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "MatLabUtilites.h"
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

double v1[] = {0.0,0.0,0.0};
double v2[]={4,-1,2};
double v3[]={2,-2,-1};


//------------------------------------------------------------------------------ 
//  int main()
//------------------------------------------------------------------------------ 
/**  
* Execute the test for MatLabUtilites.c  
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
  // Test norm
  printf("---- Test NORM ----\n");
  
  assert(fabs(norm(v1) - 0.0) < EPSILON);
  printf("  DOT función: %f\n", dot(v2,v3));
  printf("  DOT real: %f\n", 8.0);
  printf("  Diferencia: %f\n", fabs(dot(v2,v3) - 8.0));

  printf(GREEN "---- Pass Test NORM ----\n" RESET);

  //Test sign
  printf("---- Test SIGN ----\n" RESET);

  assert(fabs(sign(0) + 1.0) < EPSILON);
  printf("  DOT función: %f\n", dot(v2,v3));
  printf("  DOT real: %f\n", 8.0);
  printf("  Diferencia: %f\n", fabs(dot(v2,v3) - 8.0));

  assert(fabs(sign(-1) + 1.0) < EPSILON);
  printf("  DOT función: %f\n", dot(v2,v3));
  printf("  DOT real: %f\n", 8.0);
  printf("  Diferencia: %f\n", fabs(dot(v2,v3) - 8.0));

  assert(fabs(sign(5) - 1.0) < EPSILON);
  printf("  DOT función: %f\n", dot(v2,v3));
  printf("  DOT real: %f\n", 8.0);
  printf("  Diferencia: %f\n", fabs(dot(v2,v3) - 8.0));

  printf(GREEN "---- Pass Test SIGN ----\n" RESET);

  //Test dot
  printf("---- Test DOT ----\n" RESET);

  assert(fabs(dot(v2,v3) - 8.0) < EPSILON);
  printf("  DOT función: %f\n", dot(v2,v3));
  printf("  DOT real: %f\n", 8.0);
  printf("  Diferencia: %f\n", fabs(dot(v2,v3) - 8.0));

  printf(GREEN "---- Pass Test DOT ----\n" RESET);

  //Test traspuesta
  printf("---- Test TRASPUESTA ----\n");
  double matriz1 [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  double res [3][3];
  traspuesta(matriz1,res);
  double x[3][3]={{1,4,7},{2,5,8},{3,6,9}};
  assert(matricesIguales(res,x));
  printf(GREEN "---- Pass Test TRASPUESTA ----\n" RESET);

  //Test det
  printf("---- Test DET ----\n");

  double xDet [3][3]={{1,4,-1},{-1,3,2},{2,2,0}};
  assert(fabs(det(xDet) - 20.0) < EPSILON);

  printf("  Det función: %f\n", det(xDet));
  printf("  Det real: %f\n", 20.0);
  printf("  Diferencia: %f\n", fabs(det(xDet) - 20.0));

  printf(GREEN "---- Pass Test DET ----\n" RESET);

  //Test zeros
  printf("---- Test ZEROS ----\n");
  double zerosMatrix [3][3]={{0,0,0},{0,0,0},{0,0,0}};
  zeros(res);
  assert(matricesIguales(res,zerosMatrix));
  printf(GREEN "---- Pass Test ZEROS ----\n" RESET);

  //Test fix
  printf("---- Test FIX ----\n");
  assert(fabs(fix(-1.9) + 1.0) < EPSILON);
  assert(fabs(fix(1.6) - 1.0) < EPSILON);
  assert(fabs(fix(-4.5) + 4.0) < EPSILON);
  assert(fabs(fix(4.5) - 4.0) < EPSILON);
  printf(GREEN "---- Pass Test FIX ----\n" RESET);
  
  //Test abs
  printf("---- Test ABS ----\n");
  double positive = 3.0;
  double zero = 0.0;
  double negative = -9.0;

  assert(abs(positive) == fabs(positive));
  assert(abs(zero) == fabs(zero));
  assert(abs(negative) == fabs(negative));

  printf(GREEN "---- Pass Test ABS ----\n" RESET);

  //Test all
  printf("---- Test ALL ----\n");
  double xAll[3][3]={{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
  assert(all(xAll,-1));
  printf(GREEN "---- Pass Test ALL ----\n" RESET);

  //Test sumaMatrices
  printf("---- Test SUMA MATRICES ----\n");
  double matriz2 [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  sumaMatrices(matriz1,matriz2,res);
  double xSumaMatrices[3][3]={{2,4,6},{8,10,12},{14,16,18}};
  assert(matricesIguales(res,xSumaMatrices));
  printf(GREEN "---- Pass Test SUMA MATRICES ----\n" RESET);


  //Test restaMatrices
  printf("---- Test RESTA MATRICES ----\n");
  restaMatrices(matriz1,matriz2,res);
  double xRestaMatrices[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  assert(matricesIguales(res,xRestaMatrices));
  printf(GREEN "---- Pass Test RESTA MATRICES ----\n" RESET);

  //Test multiplicacionMatrices
  printf("---- Test MULTIPLICACION MATRICES ----\n");
  double matrixA [3][3] = {{1,-1,1},{2,2,3},{-2,-3,-1}};
  double matrixB [3][3] = {{1,0,4},{0,2,5},{1,3,0}};
  double matrixC [3][3];
  double matrixResultado [3][3] = {{2,1,-1}, {5,13,18},{-3,-9,-23}};
  multiplicacionMatrices(matrixA,matrixB,matrixC);
  assert(matricesIguales(matrixC,matrixResultado));
  printf(GREEN "---- Pass Test MULTIPLICACION MATRICES ----\n" RESET);
  
  //Test cross
  printf("---- Test CROSS ----\n");
  //Cross vectores
  double vectorRes[3];
  double vector1Cross[3]={4,-2,1};
  double vector2Cross[3]={1,-1,3};
  crossVector(vector1Cross,vector2Cross,vectorRes);
  double vectorResCross[3]={-5,-11,-2};
  assert(vectoresIguales(vectorRes,vectorResCross));

  //Cross matrices

  printf(GREEN "---- Pass Test CROSS ----\n" RESET);
  
  //Test roots
  printf("---- Test ROOTS ----\n");
  double polinome[9] = {1,0, -73120740632072, 0, 0, -1.58793679567638e+36, 0, 0, -1.19853848536909e+58};
  double resRoots[8];
  int numCoeficientes = 9;
  double solution[8] = {20488505.5958373, -16734286.9676343};
  int zeros;

  roots(polinome, numCoeficientes, resRoots, &zeros);

  printf("  Numero de soluciones reales: %d\n", zeros);
  printf("\n");

  for (int i = 0; i < zeros; i++) {
    printf("  Roots de C: %f\n", resRoots[i]);
    printf("  Roots de Matlab: %f\n", solution[i]);
    printf("  Diferencia: %f\n", fabs(resRoots[i] - solution[i]));

    assert(fabs(resRoots[i] - solution[i]) < EPSILON);

    printf("\n");
  }

  printf(GREEN "---- Pass Test ROOTS ----\n");

  //Pasa todos los test
  printf(GREEN "---- All Pass Test ----\n");
}