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

double dot(double vector1 [], double vector2 []);
double norm(double vector []);
int sign(double numero);

void traspuesta(double matriz [ROWS][COLS], double resultado [ROWS][COLS]);
bool matricesIguales(double matriz1 [ROWS][COLS], double matriz2 [ROWS][COLS]);
bool vectoresIguales(double vector1[] , double vector2[]);

double det(double matriz [ROWS][COLS]);

void zeros(double matriz [ROWS][COLS]);
void zeros_vector(double vector [ROWS]);

double fix(double num);

bool all(double matriz [ROWS][COLS], double num);
bool allVector(double vector[],double num);

void sumaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void restaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void multiplicacionMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void crossMatrix(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void crossVector(double vector1 [],double vector2 [],double resultado []);

void sumaVectores(double vector1[], double vector2[], double res[]);
void restaVectores(double vector1[], double vector2[], double res[]);
void multiplicacionVectorPorEscalar(double vector[], double valor, double res[]);
void devisionVectorPorEscalar(double vector[], double valor, double res[]);

void opuestoVector(double vector[], double res[]);

void roots(double poly[], int numCoeficientes, double solucionesReales[], int *numSolucionesReales);

void muestraVector(double vector[]);
void muestraMatriz(double matrix[ROWS][COLS]);

double matlab_mod(double a, double q);

void divideComponentesVectorEntreValor(double vector[], double valor, double res[]);
void multiplicaComponentesVectorPorValor(double vector[], double valor, double res[]);

void copiaVector(double original [], double copia[]);

#endif /* MATLABUTILITIES_H */