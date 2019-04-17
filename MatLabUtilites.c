#include <stdio.h>
#include "MatLabUtilites.h"

// Como Juan Félix dijo que los vectores de tamaño fijo para ser coherentes no es necesario poner el tamaño

double norm(double *vector) {
	return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
}

double dot(double *vector1, double *vector2) {
  return sqrt(vector1[0]*vector1[0]+vector1[1]*vector1[1]+vector1[2]*vector1[2])*sqrt(vector2[0]*vector2[0]+vector2[1]*vector2[1]+vector2[2]*vector2[2]);
}