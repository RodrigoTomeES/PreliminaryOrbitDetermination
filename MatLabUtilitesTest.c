// Test MatLabUtilites
#include <assert.h>
#include <math.h>
#include "MatLabUtilites.h"
#define EPSILON pow(10, -13)

double v1[] = {0.0,0.0,0.0};

// Test norm
int main () {
  assert(fabs(norm(v1) - 0.0) < EPSILON);
}