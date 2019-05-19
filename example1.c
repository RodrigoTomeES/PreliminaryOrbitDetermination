#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Mjday.h"
#include "Position.h"
#include "MatLabUtilites.h"
//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
* Execute the test for Preliminary Orbit Determination Proyect
* This function show the results of the tests.
*
* @param  - none
* @return - none
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
int main () {
    // read Earth orientation parameters
    FILE* fid = fopen("eop19620101.txt","rt");

    int filas = 20026;
    int columnas = 13;
    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;

    if (fid == NULL){
        exit(EXIT_FAILURE);
    }

    double **eop;
    eop = (double **) malloc (filas*sizeof(double *));

    if (eop != NULL) {
        for (int i = 0; i < filas; i++) {
            eop[i] = (double *) malloc (columnas * sizeof(double));
            if (fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) == EOF) {
                break;
            }
            eop[i][0]  =  v1;
            eop[i][1]  =  v2;
            eop[i][2]  =  v3;
            eop[i][3]  =  v4;
            eop[i][4]  =  v5;
            eop[i][5]  =  v6;
            eop[i][6]  =  v7;
            eop[i][7]  =  v8;
            eop[i][8]  =  v9;
            eop[i][9]  = v10;
            eop[i][10] = v11;
            eop[i][11] = v12;
            eop[i][12] = v13;

            //printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }
    } else {
        printf("Es null");
    }

    fclose(fid);
    
    // read observations
    fid = fopen("sat1.txt","r");

    int i = 0;

    int Y,M,D,h,m;
    double s;
    double rtasc, decl;

    double obs[3][3];

    while(i<3){
        if (fid==NULL){
            break;
        }
        
        Y=0;
        M=0;
        D=0;
        h=0;
        m=0;
        s=0.0;
        rtasc=0;
        decl=0;

        size_t len=0;
        size_t buff=0;
        char * linea=NULL;
        len=getline(&linea,&buff,fid);
        
        //Y = str2num(tline(1:4));
        char ybuff[4];
        strncpy(ybuff, linea, 4);
        Y=atoi(ybuff);
        printf("Y %d\n",Y);
        

        //M = str2num(tline(6:7));
        char mbuff[2];
        strncpy(mbuff, linea+5, 2);
        M=atoi(mbuff);
        printf("M %d\n",M);

        //D = str2num(tline(9:10));
        char dbuff[2];
        strncpy(dbuff, linea+8, 2);
        D=atoi(dbuff);
        printf("D %d\n",D);

        //h = str2num(tline(12:13));
        char hbuff[2];
        strncpy(hbuff, linea+11, 2);
        h=atoi(hbuff);
        printf("h %d\n",h);

        //m = str2num(tline(15:16));
        char mmbuff[2];
        strncpy(mmbuff, linea+14, 2);
        m=atoi(mmbuff);
        printf("mm %d\n",m);

        //s = str2num(tline(18:23));
        char sbuff[2];
        strncpy(sbuff, linea+17, 6);
        s=atof(sbuff);
        printf("s %f\n",s);


        //rtasc = str2num(tline(24:35));
        char rtascbuff[12];
        strncpy(rtascbuff, linea+23,12);
        rtasc=atof(rtascbuff);
        printf("rtasc %f\n",rtasc);



        //decl = str2num(tline(36:end));
        char declbuff[len-35];
        strncpy(declbuff, linea+35,len-35);
        decl=atof(declbuff);
        printf("decl %f\n",decl);

        //obs[i][0] = Mjday(Y,M,D,h,m,s);
        //obs[i][1] = Rad*rtasc;
        //obs[i][2] = Rad*decl;
        i = i+1;
    }

    fclose(fid);

    
}