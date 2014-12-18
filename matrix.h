//
//  matrix.h
//  svd_j
//
//  Created by hang shao on 12/13/14.
//  Copyright (c) 2014 hang shao. All rights reserved.
//

#ifndef _MATRIX_H
#define _MATRIX_H

#include<assert.h>
#include<stdlib.h>
#include<stdio.h>


double** getMatrix(int rows,int columns){
    double **rect=(double**)calloc(rows,sizeof(double*));
    int i;
    for(i=0;i<rows;++i)
        rect[i]=(double*)calloc(columns,sizeof(double));
    return rect;
}


double** getIndentityMatrix(int rows){
    double** IM=getMatrix(rows,rows);
    int i;
    for(i=0;i<rows;++i)
        IM[i][i]=1.0;
    return IM;
}


double** copyMatrix(double** matrix,int rows,int columns){
    double** rect=getMatrix(rows,columns);
    int i,j;
    for(i=0;i<rows;++i)
        for(j=0;j<columns;++j)
            rect[i][j]=matrix[i][j];
    return rect;
}


void getFromArray(double** matrix,int rows,int columns,double *arr){
    int i,j,k=0;
    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            matrix[i][j]=arr[k++];
        }
    }
}


void printMatrix(double** matrix,int rows,int columns){
    int i,j;
    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            printf("%-10f\t",matrix[i][j]);
        }
        printf("\n");
    }
}


void freeMatrix(double** matrix,int rows){
    int i;
    for(i=0;i<rows;++i)
        free(matrix[i]);
    free(matrix);
}


double* getRow(double **matrix,int rows,int columns,int index){
    assert(index<rows);
    double *rect=(double*)calloc(columns,sizeof(double));
    int i;
    for(i=0;i<columns;++i)
        rect[i]=matrix[index][i];
    return rect;
}


double* getColumn(double **matrix,int rows,int columns,int index){
    assert(index<columns);
    double *rect=(double*)calloc(rows,sizeof(double));
    int i;
    for(i=0;i<rows;++i)
        rect[i]=matrix[i][index];
    return rect;
}


void setColumn(double **matrix,int rows,int columns,int index,double *arr){
    assert(index<columns);
    int i;
    for(i=0;i<rows;++i)
        matrix[i][index]=arr[i];
}


void exchangeColumn(double **matrix,int rows,int columns,int i,int j){
    assert(i<columns);
    assert(j<columns);
    int row;
    for(row=0;row<rows;++row){
        double tmp=matrix[row][i];
        matrix[row][i]=matrix[row][j];
        matrix[row][j]=tmp;
    }
}



double** getTranspose(double **matrix,int rows,int columns){
    double **rect=getMatrix(columns,rows);
    int i,j;
    for(i=0;i<columns;++i){
        for(j=0;j<rows;++j){
            rect[i][j]=matrix[j][i];
        }
    }
    return rect;
}


double vectorProduct(double *vector1,double *vector2,int len){
    double rect=0.0;
    int i;
    for(i=0;i<len;++i)
        rect+=vector1[i]*vector2[i];
    return rect;
}


double** matrixProduct(double **matrix1,int rows1,int columns1,double **matrix2,int columns2){
    double **rect=getMatrix(rows1,columns2);
    int i,j;
    for(i=0;i<rows1;++i){
        for(j=0;j<columns2;++j){
            double *vec1=getRow(matrix1,rows1,columns1,i);
            double *vec2=getColumn(matrix2,columns1,columns2,j);
            rect[i][j]=vectorProduct(vec1,vec2,columns1);
            free(vec1);
            free(vec2);
        }
    }
    return rect;
}


double getColumnNorm(double** matrix,int rows,int columns,int index){
    assert(index<columns);
    double* vector=getColumn(matrix,rows,columns,index);
    double norm=vectorProduct(vector,vector,rows);
    free(vector);
    return norm;
}


void printVector(double* vector,int len){
    int i;
    for(i=0;i<len;++i)
        printf("%-15.8f\t",vector[i]);
    printf("\n");
}

#endif
