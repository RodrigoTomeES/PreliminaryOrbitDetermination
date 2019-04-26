#include <stdio.h>
#include "MatLabUtilites.h"

// Como Juan Félix dijo que los vectores de tamaño fijo para ser coherentes no es necesario poner el tamaño

double norm(double vector  []) {
	return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
}

double dot(double vector1 [], double vector2 []) {
  return vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2];
}

int sign(double numero){
  int res=-1;

  if(numero>0){
    res=1;
  }

  return res;
}

void traspuesta(double matriz [ROWS][COLS], double resultado [ROWS][COLS]){
  resultado[0][0]=matriz[0][0];
  resultado[0][1]=matriz[1][0];
  resultado[0][2]=matriz[2][0];
  resultado[1][0]=matriz[0][1];
  resultado[1][1]=matriz[1][1];
  resultado[1][2]=matriz[2][1];
  resultado[2][0]=matriz[0][2];
  resultado[2][1]=matriz[1][2];
  resultado[2][2]=matriz[2][2];
}

bool matricesIguales(double matriz1 [ROWS][COLS], double matriz2 [ROWS][COLS]){
  return (matriz1[0][0]==matriz2[0][0]&&matriz1[0][1]==matriz2[0][1]&&matriz1[0][2]==matriz2[0][2]&&matriz1[1][0]==matriz2[1][0]&&matriz1[1][1]==matriz2[1][1]&&matriz1[1][2]==matriz2[1][2]&&matriz1[2][0]==matriz2[2][0]&&matriz1[2][1]==matriz2[2][1]&&matriz1[2][2]==matriz2[2][2]);
}

double det(double matriz [ROWS][COLS]){
  return matriz[0][0]*matriz[1][1]*matriz[2][2]+matriz[0][2]*matriz[1][0]*matriz[2][1]+matriz[0][1]*matriz[1][2]*matriz[2][0]-(matriz[0][2]*matriz[1][1]*matriz[2][0]+matriz[1][0]*matriz[0][1]*matriz[2][2]+matriz[0][0]*matriz[1][2]*matriz[2][1]);
}

void zeros(double matriz [ROWS][COLS]){
  matriz[0][0]=0;
  matriz[0][1]=0;
  matriz[0][2]=0;
  matriz[1][0]=0;
  matriz[1][1]=0;
  matriz[1][2]=0;
  matriz[2][0]=0;
  matriz[2][1]=0;
  matriz[2][2]=0;
}

double fix(double num){
  double res = 0;

  if(sign(num)==1){
    res=floor(num);
  }else{
    res=ceil(num);
  }

  return res;
}
double abs(double num){
  return fabs(num);
}

bool all(double matriz [ROWS][COLS], double num){
  return matriz[0][0]==num&&
          matriz[0][1]==num&&
          matriz[0][2]==num&&
          matriz[1][0]==num&&
          matriz[1][1]==num&&
          matriz[1][2]==num&&
          matriz[2][0]==num&&
          matriz[2][1]==num&&
          matriz[2][2]==num;
}

void sumaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]){
  resultado[0][0]=matriz1[0][0]+matriz2[0][0];
  resultado[0][1]=matriz1[0][1]+matriz2[0][1];
  resultado[0][2]=matriz1[0][2]+matriz2[0][2];
  resultado[1][0]=matriz1[1][0]+matriz2[1][0];
  resultado[1][1]=matriz1[1][1]+matriz2[1][1];
  resultado[1][2]=matriz1[1][2]+matriz2[1][2];
  resultado[2][0]=matriz1[2][0]+matriz2[2][0];
  resultado[2][1]=matriz1[2][1]+matriz2[2][1];
  resultado[2][2]=matriz1[2][2]+matriz2[2][2];
}

void restaMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]){
  resultado[0][0]=matriz1[0][0]-matriz2[0][0];
  resultado[0][1]=matriz1[0][1]-matriz2[0][1];
  resultado[0][2]=matriz1[0][2]-matriz2[0][2];
  resultado[1][0]=matriz1[1][0]-matriz2[1][0];
  resultado[1][1]=matriz1[1][1]-matriz2[1][1];
  resultado[1][2]=matriz1[1][2]-matriz2[1][2];
  resultado[2][0]=matriz1[2][0]-matriz2[2][0];
  resultado[2][1]=matriz1[2][1]-matriz2[2][1];
  resultado[2][2]=matriz1[2][2]-matriz2[2][2];
}

void multiplicacionMatrices(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]){

}

void crossMatrix(double matriz1 [ROWS][COLS],double matriz2 [ROWS][COLS],double resultado [ROWS][COLS]){

}

bool vectoresIguales(double vector1[] , double vector2[]){
  return vector1[0]==vector2[0]&&vector1[1]==vector2[1]&&vector1[2]==vector2[2];
}

void crossVector(double matriz1 [],double matriz2 [],double resultado []){
  resultado[0]=matriz1[1]*matriz2[2]-matriz1[2]*matriz2[1];
  resultado[1]=matriz1[2]*matriz2[0]-matriz1[0]*matriz2[2];
  resultado[2]=matriz1[0]*matriz2[1]-matriz1[1]*matriz2[0];
}

void roots(double vector [], int num, double res[]){

}