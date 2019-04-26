#include <math.h>
#include <stdbool.h>

#define ROWS 3
#define COLS 3

double dot(double vector1 [], double vector2 []);
double norm(double vector []);
int sign(double numero);

void traspuesta(double matriz [ROWS][COLS], double resultado [ROWS][COLS]);
bool matricesIguales(double matriz1 [ROWS][COLS], double matriz2 [ROWS][COLS]);

double det(double matriz [ROWS][COLS]);

void zeros(double matriz [ROWS][COLS]);

double fix(double num);
double abs(double num);

bool all(double matriz [ROWS][COLS], double num);

void sumaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void restaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void multiplicacionMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]);
void cross(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS]);

void roots(double poly[], int numCoeficientes, double solucionesReales[], int *numSolucionesReales);