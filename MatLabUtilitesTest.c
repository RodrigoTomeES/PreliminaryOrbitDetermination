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
  double matriz1 [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  double res [3][3];
  traspuesta(matriz1,res);
  double x[3][3]={{1,4,7},{2,5,8},{3,6,9}};
  assert(matricesIguales(res,x));

  //Test det
  
  double xDet [3][3]={{1,4,-1},{-1,3,2},{2,2,0}};
  assert(fabs(det(xDet) - 20.0) < EPSILON);

  //Test zeros
  
  //Test fix
  
  //Test abs

  //Test all
  double xAll[3][3]={{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
  assert(all(xAll,-1));

  //Test sumaMatrices
  double matriz2 [3][3]={{1,2,3},{4,5,6},{7,8,9}};
  sumaMatrices(matriz1,matriz2,res);
  double xSumaMatrices[3][3]={{2,4,6},{8,10,12},{14,16,18}};
  assert(matricesIguales(res,xSumaMatrices));


  //Test restaMatrices
  restaMatrices(matriz1,matriz2,res);
  double xRestaMatrices[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  assert(matricesIguales(res,xRestaMatrices));

  //Test multiplicacionMatrices
  
  //Test cross
  
  //Test roots

}