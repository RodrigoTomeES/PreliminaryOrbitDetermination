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

  //Test sign
  assert(fabs(sign(0) + 1.0) < EPSILON);
  assert(fabs(sign(-1) + 1.0) < EPSILON);
  assert(fabs(sign(5) - 1.0) < EPSILON);

  //Test dot
  assert(fabs(dot(v2,v3) - 8.0) < EPSILON);

  //Test traspuesta
  double matriz [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  double res [3][3];
  traspuesta(matriz,res);
  double x[3][3]={{1,4,7},{2,5,8},{3,6,9}};
  assert(matricesIguales(res,x));

  //Test det
  
  //Test zeros
  double zerosMatrix [3][3]={{0,0,0},{0,0,0},{0,0,0}};
  zeros(res);
  assert(matricesIguales(res,zerosMatrix));
  printf("---- Pass Test ZEROS ----");

  //Test fix
  
  //Test abs
  double positive = 3.0;
  double zero = 0.0;
  double negative = -9.0;

  assert(abs(positive) == fabs(positive));
  assert(abs(zero) == fabs(zero));
  assert(abs(negative) == fabs(negative));

  printf("---- Pass Test ABS ----");
  //Test all
  
  //Test sumaMatrices
  
  //Test restaMatrices
  
  //Test multiplicacionMatrices
  
  //Test cross
  
  //Test roots
  
}