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
* Provides the definitions of Matlab functions that are not in C.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef MATLABUTILITIES_H
#define MATLABUTILITIES_H

#include <math.h>
#include <stdbool.h>
#include "SAT_Const.h"

#define ROWS 3
#define COLS 3

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
double dot(double vector1 [], double vector2 []);

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
double norm(double vector []);

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
int sign(double numero);

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
void traspuesta(double matriz [ROWS][COLS], double resultado [ROWS][COLS]);

//------------------------------------------------------------------------------
//  void traspuesta_vector(double vector [], double resultado [1][COLS])
//------------------------------------------------------------------------------
/**
* Performs the transpose of a vector
*
* @param  - double vector [], double resultado [1][COLS]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void traspuesta_vector(double vector [], double resultado [1][COLS]);

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
bool matricesIguales(double matriz1 [ROWS][COLS], double matriz2 [ROWS][COLS]);

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
bool vectoresIguales(double vector1[] , double vector2[]);

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
double det(double matriz [ROWS][COLS]);

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
void zeros(double matriz [ROWS][COLS]);

//------------------------------------------------------------------------------
//  void zeros_vector(double vector[ROWS])
//------------------------------------------------------------------------------
/**
* Given a vector it fills its zeros components
*
* @param  - double vector[ROWS]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void zeros_vector(double vector [ROWS]);

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
double fix(double num);

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
bool all(double matriz [ROWS][COLS], double num);
bool allVector(double vector[],double num);

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
void sumaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);

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
void restaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);

//------------------------------------------------------------------------------
//  void crossMatrix(double matriz1[ROWS][COLS],
//                    double matriz2[ROWS][COLS], double resultado[ROWS][COLS])
//------------------------------------------------------------------------------
/**
* Make the cross of two matrix
*
* @param  - double matriz1[ROWS][COLS], double matriz2[ROWS][COLS], double resultado[ROWS][COLS]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void crossMatrix(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);

//------------------------------------------------------------------------------
//  void crossMatrixVector(double matriz [ROWS][COLS],double vector [],
//                         double resultado [])
//------------------------------------------------------------------------------
/**
* Make the cross of two matrix
*
* @param  - double matriz [ROWS][COLS],double vector [],double resultado []
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void crossMatrixVector(double matriz [ROWS][COLS],double vector [],double resultado []);

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
void crossVector(double vector1 [],double vector2 [],double resultado []);

//------------------------------------------------------------------------------
//  void sumaVectores(double vector1[], double vector2[], double res[])
//------------------------------------------------------------------------------
/**
* Make the plus of two vectors
*
* @param  - double vector1[], double vector2[], double res[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void sumaVectores(double vector1[], double vector2[], double res[]);

//------------------------------------------------------------------------------
//  void restaVectores(double vector1[], double vector2[], double res[])
//------------------------------------------------------------------------------
/**
* Make the minus of two vectors
*
* @param  - double vector1[], double vector2[], double res[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void restaVectores(double vector1[], double vector2[], double res[]);

//------------------------------------------------------------------------------
//  void multiplicacionVectorPorEscalar(double vector[], double valor, double res[])
//------------------------------------------------------------------------------
/**
* Make the product of a vector to climb
*
* @param  - double vector[], double valor, double res[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void multiplicacionVectorPorEscalar(double vector[], double valor, double res[]);

//------------------------------------------------------------------------------
//  void divisionVectorPorEscalar(double vector[], double valor, double res[])
//------------------------------------------------------------------------------
/**
* Make the divison of a vector to climb
*
* @param  - double vector[], double valor, double res[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void divisionVectorPorEscalar(double vector[], double valor, double res[]);

//------------------------------------------------------------------------------
//  void opuestoVector(double vector[], double res[])
//------------------------------------------------------------------------------
/**
* Does the opposite of a vector
*
* @param  - double vector[], double res[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void opuestoVector(double vector[], double res[]);

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
void roots(double poly[], int numCoeficientes, double solucionesReales[], int *numSolucionesReales);

//------------------------------------------------------------------------------
//  void muestraVector(double vector[])
//------------------------------------------------------------------------------
/**
* Show the vector
*
* @param  - double vector[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void muestraVector(double vector[]);

//------------------------------------------------------------------------------
//  void muestraMatriz(double matrix[ROWS][COLS])
//------------------------------------------------------------------------------
/**
* Show the maxtrix
*
* @param  - double matrix[ROWS][COLS]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void muestraMatriz(double matrix[ROWS][COLS]);

//------------------------------------------------------------------------------
//  double matlab_mod(double a, double q)
//------------------------------------------------------------------------------
/**
* Show the maxtrix
*
* @param  - double a, double q
* @return - double
* @exception - none
* @see - none
* @note - Mod in C doesn't equal than Matlab
*         https://stackoverflow.com/questions/28888619/modulo-function-in-c-that-behaves-like-mod-in-matlab
*/
//------------------------------------------------------------------------------
double matlab_mod(double a, double q);

//------------------------------------------------------------------------------
//  void copiaVector(double original [], double copia[])
//------------------------------------------------------------------------------
/**
* Copy the vector original in copy
*
* @param  - double original [], double copia[]
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void copiaVector(double original [], double copia[]);

#endif /* MATLABUTILITIES_H */