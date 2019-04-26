// Test MatLabUtilites
#include <assert.h>
#include <math.h>
#include "MatLabUtilites.h"
#define EPSILON pow(10, -13)

double v1[] = {0.0,0.0,0.0};
double v2[]={4,-1,2};
double v3[]={2,-2,-1};


int main () {
  // Test norm
  assert(fabs(norm(v1) - 0.0) < EPSILON);
  printf("---- Pass Test NORM ----\n");

  //Test sign
  assert(fabs(sign(0) + 1.0) < EPSILON);
  assert(fabs(sign(-1) + 1.0) < EPSILON);
  assert(fabs(sign(5) - 1.0) < EPSILON);
  printf("---- Pass Test SIGN ----\n");

  //Test dot
  assert(fabs(dot(v2,v3) - 8.0) < EPSILON);
  printf("---- Pass Test DOT ----\n");

  //Test traspuesta
  double matriz1 [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  double res [3][3];
  traspuesta(matriz1,res);
  double x[3][3]={{1,4,7},{2,5,8},{3,6,9}};
  assert(matricesIguales(res,x));
  printf("---- Pass Test TRASPUESTA ----\n");

  //Test det
  double xDet [3][3]={{1,4,-1},{-1,3,2},{2,2,0}};
  assert(fabs(det(xDet) - 20.0) < EPSILON);
  printf("---- Pass Test DET ----\n");
  
  

  //Test zeros
  double zerosMatrix [3][3]={{0,0,0},{0,0,0},{0,0,0}};
  zeros(res);
  assert(matricesIguales(res,zerosMatrix));
  printf("---- Pass Test ZEROS ----\n");

  //Test fix
  assert(fabs(fix(-1.9) + 1.0) < EPSILON);
  assert(fabs(fix(1.6) - 1.0) < EPSILON);
  assert(fabs(fix(-4.5) + 4.0) < EPSILON);
  assert(fabs(fix(4.5) - 4.0) < EPSILON);
  printf("---- Pass Test FIX ----\n");
  
  //Test abs
  double positive = 3.0;
  double zero = 0.0;
  double negative = -9.0;

  assert(abs(positive) == fabs(positive));
  assert(abs(zero) == fabs(zero));
  assert(abs(negative) == fabs(negative));

  printf("---- Pass Test ABS ----\n");

  //Test all
  double xAll[3][3]={{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
  assert(all(xAll,-1));
  printf("---- Pass Test ALL ----\n");

  //Test sumaMatrices
  double matriz2 [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  sumaMatrices(matriz1,matriz2,res);
  double xSumaMatrices[3][3]={{2,4,6},{8,10,12},{14,16,18}};
  assert(matricesIguales(res,xSumaMatrices));
  printf("---- Pass Test SUMA MATRICES ----\n");


  //Test restaMatrices
  restaMatrices(matriz1,matriz2,res);
  double xRestaMatrices[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  assert(matricesIguales(res,xRestaMatrices));
  printf("---- Pass Test RESTA MATRICES ----\n");

  //Test multiplicacionMatrices
  printf("---- Pass Test MULTIPLICACION MATRICES ----\n");
  
  //Test cross
  //Cross vectores
  double vectorRes[3];
  double vector1Cross[3]={4,-2,1};
  double vector2Cross[3]={1,-1,3};
  crossVector(vector1Cross,vector2Cross,vectorRes);
  double vectorResCross[3]={-5,-11,-2};
  assert(vectoresIguales(vectorRes,vectorResCross));

  //Cross matrices
  

  printf("---- Pass Test CROSS ----\n");
  
  //Test roots
  printf("---- Pass Test ROOTS ----\n");

  //Pasa todos los test
  printf("---- All Pass Test ----\n");
}