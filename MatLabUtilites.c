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
* Provides a basic implementation of Matlab functions that are not in C.  
*  
* @note     
*/
//------------------------------------------------------------------------------

#include <stdio.h>
#include "MatLabUtilites.h"
#include "rpoly.h"
#define EPSILON pow(10, -7)

// Como Juan Félix dijo que los vectores son de tamaño fijo para ser coherentes no es necesario poner el tamaño

//------------------------------------------------------------------------------
//  double norm(double vector[])
//------------------------------------------------------------------------------
/**  
* Execute the test for MatLabUtilites.c  
* This function show the results of the tests.
*  
* @param  - double vector[]
* @return - double
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
double norm(double vector[])
{
  return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

//------------------------------------------------------------------------------
//  double dot(double vector1[], double vector2[])
//------------------------------------------------------------------------------
/**  
* Execute the test for MatLabUtilites.c  
* This function show the results of the tests.
*  
* @param  - double vector1[], double vector2[]
* @return - double
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
double dot(double vector1[], double vector2[])
{
  return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}

//------------------------------------------------------------------------------
//  int sign(double numero)
//------------------------------------------------------------------------------
/**  
* Evaluate the sign of the number and return a -1 when it is negative 
* (or zero) and a 1 when it is positive.
*  
* @param  - double numero
* @return - int
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
int sign(double numero)
{
  int res = -1;

  if (numero > 0)
  {
    res = 1;
  }

  return res;
}

//------------------------------------------------------------------------------
//  void traspuesta(double matriz[ROWS][COLS], double resultado[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Performs the transpose of a matrix
*  
* @param  - double matriz[ROWS][COLS], double resultado[ROWS][COLS]
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void traspuesta(double matriz[ROWS][COLS], double resultado[ROWS][COLS])
{
  resultado[0][0] = matriz[0][0];
  resultado[0][1] = matriz[1][0];
  resultado[0][2] = matriz[2][0];
  resultado[1][0] = matriz[0][1];
  resultado[1][1] = matriz[1][1];
  resultado[1][2] = matriz[2][1];
  resultado[2][0] = matriz[0][2];
  resultado[2][1] = matriz[1][2];
  resultado[2][2] = matriz[2][2];
}

//------------------------------------------------------------------------------
//  bool matricesIguales(double matriz1[ROWS][COLS], double matriz2[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Given two matrices it says if they are equal component to component
*  
* @param  - double matriz1[ROWS][COLS], double matriz2[ROWS][COLS]
* @return - bool
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
bool matricesIguales(double matriz1[ROWS][COLS], double matriz2[ROWS][COLS])
{
  return (matriz1[0][0] == matriz2[0][0] && matriz1[0][1] == matriz2[0][1] && matriz1[0][2] == matriz2[0][2] && matriz1[1][0] == matriz2[1][0] && matriz1[1][1] == matriz2[1][1] && matriz1[1][2] == matriz2[1][2] && matriz1[2][0] == matriz2[2][0] && matriz1[2][1] == matriz2[2][1] && matriz1[2][2] == matriz2[2][2]);
}

//------------------------------------------------------------------------------
//  double det(double matriz[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Performs the determinant of the matrix passed by parameters
*  
* @param  - double matriz[ROWS][COLS]
* @return - double
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
double det(double matriz[ROWS][COLS])
{
  return matriz[0][0] * matriz[1][1] * matriz[2][2] + matriz[0][2] * matriz[1][0] * matriz[2][1] + matriz[0][1] * matriz[1][2] * matriz[2][0] - (matriz[0][2] * matriz[1][1] * matriz[2][0] + matriz[1][0] * matriz[0][1] * matriz[2][2] + matriz[0][0] * matriz[1][2] * matriz[2][1]);
}

//------------------------------------------------------------------------------
//  void zeros(double matriz[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Given a matrix it fills its zeros components
*  
* @param  - double matriz[ROWS][COLS]
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void zeros(double matriz[ROWS][COLS])
{
  matriz[0][0] = 0;
  matriz[0][1] = 0;
  matriz[0][2] = 0;
  matriz[1][0] = 0;
  matriz[1][1] = 0;
  matriz[1][2] = 0;
  matriz[2][0] = 0;
  matriz[2][1] = 0;
  matriz[2][2] = 0;
}

//------------------------------------------------------------------------------
//  double fix(double num)
//------------------------------------------------------------------------------
/**  
* fix(X) rounds each element of X to the nearest integer toward zero. 
* For positive X, the behavior of fix is the same as floor. 
* For negative X, the behavior of fix is the same as ceil.
*  
* @param  - double num
* @return - double
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
double fix(double num)
{
  double res = 0;

  if (sign(num) == 1)
  {
    res = floor(num);
  }
  else
  {
    res = ceil(num);
  }

  return res;
}

//------------------------------------------------------------------------------
//  double abs(double num)
//------------------------------------------------------------------------------
/**  
* Make the abs of the numer
*  
* @param  - double num
* @return - double
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
double abs(double num)
{
  return fabs(num);
}

//------------------------------------------------------------------------------
//  bool all(double matriz[ROWS][COLS], double num)
//------------------------------------------------------------------------------
/**  
* Check if all the components of a matrix are equal to the number passed 
* by parameter
*  
* @param  - double matriz[ROWS][COLS], double num
* @return - bool
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
bool all(double matriz[ROWS][COLS], double num)
{
  return matriz[0][0] == num &&
         matriz[0][1] == num &&
         matriz[0][2] == num &&
         matriz[1][0] == num &&
         matriz[1][1] == num &&
         matriz[1][2] == num &&
         matriz[2][0] == num &&
         matriz[2][1] == num &&
         matriz[2][2] == num;
}

//------------------------------------------------------------------------------
//  void sumaMatrices(double matriz1[ROWS][COLS], 
//                    double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Make the sum of two matrix
*  
* @param  - double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS]
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void sumaMatrices(double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
{
  resultado[0][0] = matriz1[0][0] + matriz2[0][0];
  resultado[0][1] = matriz1[0][1] + matriz2[0][1];
  resultado[0][2] = matriz1[0][2] + matriz2[0][2];
  resultado[1][0] = matriz1[1][0] + matriz2[1][0];
  resultado[1][1] = matriz1[1][1] + matriz2[1][1];
  resultado[1][2] = matriz1[1][2] + matriz2[1][2];
  resultado[2][0] = matriz1[2][0] + matriz2[2][0];
  resultado[2][1] = matriz1[2][1] + matriz2[2][1];
  resultado[2][2] = matriz1[2][2] + matriz2[2][2];
}

//------------------------------------------------------------------------------
//  void restaMatrices(double matriz1[ROWS][COLS], 
//                     double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Make the subtraction of two matrix
*  
* @param  - double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS]
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void restaMatrices(double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
{
  resultado[0][0] = matriz1[0][0] - matriz2[0][0];
  resultado[0][1] = matriz1[0][1] - matriz2[0][1];
  resultado[0][2] = matriz1[0][2] - matriz2[0][2];
  resultado[1][0] = matriz1[1][0] - matriz2[1][0];
  resultado[1][1] = matriz1[1][1] - matriz2[1][1];
  resultado[1][2] = matriz1[1][2] - matriz2[1][2];
  resultado[2][0] = matriz1[2][0] - matriz2[2][0];
  resultado[2][1] = matriz1[2][1] - matriz2[2][1];
  resultado[2][2] = matriz1[2][2] - matriz2[2][2];
}

//------------------------------------------------------------------------------
//  void multiplicacionMatrices(double matriz1[ROWS][COLS], 
//                    double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
//------------------------------------------------------------------------------
/**  
* Make the multiplication of two matrix
*  
* @param  - double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS]
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void multiplicacionMatrices(double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
{
  resultado[0][0] = matriz1[0][0] * matriz2[0][0] + matriz1[0][1] * matriz2[1][0] + matriz1[0][2] * matriz2[2][0];
  resultado[0][1] = matriz1[0][0] * matriz2[0][1] + matriz1[0][1] * matriz2[1][1] + matriz1[0][2] * matriz2[2][1];
  resultado[0][2] = matriz1[0][0] * matriz2[0][2] + matriz1[0][1] * matriz2[1][2] + matriz1[0][2] * matriz2[2][2];

  resultado[1][0] = matriz1[1][0] * matriz2[0][0] + matriz1[1][1] * matriz2[1][0] + matriz1[1][2] * matriz2[2][0];
  resultado[1][1] = matriz1[1][0] * matriz2[0][1] + matriz1[1][1] * matriz2[1][1] + matriz1[1][2] * matriz2[2][1];
  resultado[1][2] = matriz1[1][0] * matriz2[0][2] + matriz1[1][1] * matriz2[1][2] + matriz1[1][2] * matriz2[2][2];

  resultado[2][0] = matriz1[2][0] * matriz2[0][0] + matriz1[2][1] * matriz2[1][0] + matriz1[2][2] * matriz2[2][0];
  resultado[2][1] = matriz1[2][0] * matriz2[0][1] + matriz1[2][1] * matriz2[1][1] + matriz1[2][2] * matriz2[2][1];
  resultado[2][2] = matriz1[2][0] * matriz2[0][2] + matriz1[2][1] * matriz2[1][2] + matriz1[2][2] * matriz2[2][2];
}

// void crossMatrix(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]){

// }

//------------------------------------------------------------------------------
//  void roots(double poly[], int numCoeficientes, double solucionesReales[], int *numSolucionesReales)
//------------------------------------------------------------------------------
/**  
* Find the real roots of a polynomial of degree N
*  
* @param  - double poly[], int numCoeficientes, double solucionesReales[], int *numSolucionesReales
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void roots(double poly[], int numCoeficientes, double solucionesReales[], int *numSolucionesReales)
{
  int degree = numCoeficientes - 1;
  double real[degree];
  double imaginary[degree];
  int zeros = real_poly_roots(poly, degree, real, imaginary);
  int contadorSolucionesReales = 0;

  for (int i = 0; i < zeros; i++)
  {
    if (fabs(imaginary[i]) <= EPSILON || imaginary[i] == 0.0)
    {
      solucionesReales[i] = real[i];
      contadorSolucionesReales += 1;
    }
  }
  *numSolucionesReales = contadorSolucionesReales;
}

//------------------------------------------------------------------------------
//  bool vectoresIguales(double vector1[], double vector2[])
//------------------------------------------------------------------------------
/**  
* Compare two component-to-component vectors
*  
* @param  - double vector1[], double vector2[]
* @return - bool
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
bool vectoresIguales(double vector1[], double vector2[])
{
  return vector1[0] == vector2[0] && vector1[1] == vector2[1] && vector1[2] == vector2[2];
}

//------------------------------------------------------------------------------
//  void crossVector(double matriz1[], double matriz2[], double resultado[])
//------------------------------------------------------------------------------
/**  
* Make the cross product of two matrices
*  
* @param  - double matriz1[], double matriz2[], double resultado[]
* @return - none
* @exception - none
* @see - none  
* @note - none  
*/
//------------------------------------------------------------------------------
void crossVector(double matriz1[], double matriz2[], double resultado[])
{
  resultado[0] = matriz1[1] * matriz2[2] - matriz1[2] * matriz2[1];
  resultado[1] = matriz1[2] * matriz2[0] - matriz1[0] * matriz2[2];
  resultado[2] = matriz1[0] * matriz2[1] - matriz1[1] * matriz2[0];
}

void sumaVectores(double vector1[], double vector2[], double res[]){
  res[0]=vector1 [0]+vector2[0];
  res[1]=vector1 [1]+vector2[1];
  res[2]=vector1 [2]+vector2[2];
}

void multiplicacionVectorPorEscalar(double vector[], double valor, double res[]){
  res[0]=vector[0]*valor;
  res[1]=vector[1]*valor;
  res[2]=vector[2]*valor;
}

void devisionVectorPorEscalar(double vector[], double valor, double res[]){
  res[0]=vector[0]/valor;
  res[1]=vector[1]/valor;
  res[2]=vector[2]/valor;
}

void opuestoVector(double vector[], double res[]){
  res[0]=-vector[0];
  res[1]=-vector[1];
  res[2]=-vector[2];
}