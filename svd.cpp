//
//  main.cpp
//  svd_j
//
//  Created by hang shao on 12/13/14.
//  Copyright (c) 2014 hang shao. All rights reserved.
//

#include"mpi.h"
#include"matrix.h"
#include<string.h>
#include<stdlib.h>
#include<math.h>      //mpic++ compile with -lm

#define THREASHOLD 1e-8
#define ITERATION 20
#define  ROW 3
#define COL 3
#define LINELEN 5*COL       //length of each line for input matrix file

int sign(double number) {
    if(number<0)
        return -1;
    else
        return 1;
}

int myRank;
int procSize;
MPI_Status status;


void readFromFile(double **matrix,int row,int col,char* file){
    FILE *fp;
    int len=col*10;
    char *buf=(char*)calloc(len,sizeof(char));
    if((fp=fopen(file,"r"))==NULL){
        perror("fopen");
        printf("%s\n",file);
        exit(1);
    }
    int i,j;
    for(i=0;i<row;++i){
        if(fgets(buf,len,fp)==NULL){
            fprintf(stderr,"文件的行数小于矩阵需要的行数\n");
            exit(1);
        }
        char *seg=strtok(buf,"\t");
        double ele=atof(seg);
        matrix[i][0]=ele;
        for(j=1;j<col;++j){
            if((seg=strtok(NULL,"\t"))==NULL){
                fprintf(stderr,"文件的列数小于矩阵需要的列数\n");
                exit(1);
            }
            ele=atof(seg);
            matrix[i][j]=ele;
        }
        memset(buf,0x00,len);
    }
    free(buf);
    fclose(fp);
}


void writeToFile(double **matrix,int rows,int columns,char* file){
    FILE *fp;
    if((fp=fopen(file,"w"))==NULL){
        perror("fopen");
        exit(1);
    }
    fprintf(fp,"%d\t%d\n",rows,columns);
    int i,j;
    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            fprintf(fp,"%-10f\t",matrix[i][j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


void vectorToFile(double *vector,int len,char* file){
    FILE *fp;
    if((fp=fopen(file,"w"))==NULL){
        perror("fopen");
        exit(1);
    }
    int i;
    for(i=0;i<len;++i){
        fprintf(fp,"%-10f\t",vector[i]);
    }
    fclose(fp);
}


void orthogonalVector(double *Ci,double *Cj,int len1,double *Vi,double *Vj,int len2,int *pass){
    double ele=vectorProduct(Ci,Cj,len1);
    if(fabs(ele)<THREASHOLD)
        return;
    *pass=0;
    double ele1=vectorProduct(Ci,Ci,len1);
    double ele2=vectorProduct(Cj,Cj,len1);
    
    
    
    double tao=(ele1-ele2)/(2*ele);
    double tan=sign(tao)/(fabs(tao)+sqrt(1+pow(tao,2)));
    double cos=1/sqrt(1+pow(tan,2));
    double sin=cos*tan;
    
    int row;
    for(row=0;row<len1;++row){
        double var1=Ci[row]*cos+Cj[row]*sin;
        double var2=Cj[row]*cos-Ci[row]*sin;
        
        Ci[row]=var1;
        Cj[row]=var2;
        
    }
    for(row=0;row<len2;++row){
        double var1=Vi[row]*cos+Vj[row]*sin;
        double var2=Vj[row]*cos-Vi[row]*sin;
        
        Vi[row]=var1;
        Vj[row]=var2;
        
    }
}


void orthogonal(double **matrix,int rows,int columns,int i,int j,int *pass,double **V){
    assert(i<j);
    
    double* Ci=getColumn(matrix,rows,columns,i);
    double* Cj=getColumn(matrix,rows,columns,j);
    double* Vi=getColumn(V,columns,columns,i);
    double* Vj=getColumn(V,columns,columns,j);
    orthogonalVector(Ci,Cj,rows,Vi,Vj,columns,pass);
    
    int row;
    for(row=0;row<rows;++row){
        matrix[row][i]=Ci[row];
        matrix[row][j]=Cj[row];
    }
    for(row=0;row<columns;++row){
        V[row][i]=Vi[row];
        V[row][j]=Vj[row];
    }
    free(Ci);
    free(Cj);
    free(Vi);
    free(Vj);
}

void normalize(double **A,int rows,int columns){
    double *sigular=(double*)calloc(columns,sizeof(double));
    int i,j;
    for(i=0;i<columns;++i){
        double *vector=getColumn(A,rows,columns,i);
        double norm=sqrt(vectorProduct(vector,vector,rows));
        sigular[i]=norm;
    }
    char outFileS[7]={'S','X','.','m','a','t','\0'};
    outFileS[1]='0'+myRank;
    vectorToFile(sigular,columns,outFileS);
    double **U=getMatrix(rows,columns);
    for(j=0;j<columns;++j){
        if(sigular[j]==0)
            for(i=0;i<rows;++i)
                U[i][j]=0;
        else
            for(i=0;i<rows;++i)
                U[i][j]=A[i][j]/sigular[j];
    }
    char outFileU[7]={'U','X','.','m','a','t','\0'};
    outFileU[1]='0'+myRank;
    writeToFile(U,rows,columns,outFileU);
    free(sigular);
    freeMatrix(U,rows);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&procSize);
    
    assert(myRank<10);
    int totalColumn=COL*procSize;
    
    /*准备矩阵块A和V*/
    char matrixFile[11]={'b','l','o','c','k','X','.','m','a','t','\0'};
    matrixFile[5]='0'+myRank;
    double **A=getMatrix(ROW,COL);
    readFromFile(A,ROW,COL,matrixFile);
    double **V=getMatrix(totalColumn,COL);
    int j;
    for(j=0;j<COL;++j){
        V[COL*myRank+j][j]=1.0;
    }
    
    
    int iteration=ITERATION;
    while(iteration-->0){
        
        int pass=1;
        int allpass=0;
        int i;
        for(i=1;i<=totalColumn;++i){
            int j;
            int send=0,recv=0;
            int mod1=i%2;
            int mod2=(myRank*COL)%2;
            if(mod1^mod2){
                j=0;
            }
            else{
                j=1;
                if(myRank>0){
                    recv=1;
                }
            }
            for(++j;j<COL;j+=2){
                orthogonal(A,ROW,COL,j-1,j,&pass,V);
                exchangeColumn(A,ROW,COL,j-1,j);
                exchangeColumn(V,totalColumn,COL,j-1,j);
            }
            if(j==COL){
                if(myRank<procSize-1){
                    send=1;
                }
            }
            if(send){
                double *lastColumnA=getColumn(A,ROW,COL,COL-1);
                double *lastColumnV=getColumn(V,totalColumn,COL,COL-1);
                MPI_Send(lastColumnA,ROW,MPI_DOUBLE,myRank+1,59,MPI_COMM_WORLD);
                MPI_Send(lastColumnV,totalColumn,MPI_DOUBLE,myRank+1,60,MPI_COMM_WORLD);
                free(lastColumnA);
                free(lastColumnV);
            }
            if(recv){
                double* preColumnA=(double*)calloc(ROW,sizeof(double));
                double* preColumnV=(double*)calloc(totalColumn,sizeof(double));
                MPI_Recv(preColumnA,ROW,MPI_DOUBLE,myRank-1,59,MPI_COMM_WORLD,&status);
                MPI_Recv(preColumnV,totalColumn,MPI_DOUBLE,myRank-1,60,MPI_COMM_WORLD,&status);
                double* firstColumnA=getColumn(A,ROW,COL,0);
                double* firstColumnV=getColumn(V,totalColumn,COL,0);
                orthogonalVector(preColumnA,firstColumnA,ROW,preColumnV,firstColumnV,totalColumn,&pass);
                setColumn(A,ROW,COL,0,preColumnA);
                setColumn(V,totalColumn,COL,0,preColumnV);
                MPI_Send(firstColumnA,ROW,MPI_DOUBLE,myRank-1,49,MPI_COMM_WORLD);
                MPI_Send(firstColumnV,totalColumn,MPI_DOUBLE,myRank-1,50,MPI_COMM_WORLD);
                free(preColumnA);
                free(preColumnV);
                free(firstColumnA);
                free(firstColumnV);
            }
            if(send){
                double* nextColumnA=(double*)calloc(ROW,sizeof(double));
                double* nextColumnV=(double*)calloc(totalColumn,sizeof(double));
                MPI_Recv(nextColumnA,ROW,MPI_DOUBLE,myRank+1,49,MPI_COMM_WORLD,&status);
                MPI_Recv(nextColumnV,totalColumn,MPI_DOUBLE,myRank+1,50,MPI_COMM_WORLD,&status);
                setColumn(A,ROW,COL,COL-1,nextColumnA);
                setColumn(V,totalColumn,COL,COL-1,nextColumnV);
                free(nextColumnA);
                free(nextColumnV);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&pass,&allpass,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Bcast(&allpass,1,MPI_INT,0,MPI_COMM_WORLD);
        if(allpass==procSize)
            break;
    }
    if(myRank==0){
        printf("迭代次数:%d\n",ITERATION-iteration-1);
    }
    char outFileV[7]={'V','X','.','m','a','t','\0'};
    outFileV[1]='0'+myRank;
    writeToFile(V,totalColumn,COL,outFileV);
    normalize(A,ROW,COL);
    freeMatrix(A,ROW);
    freeMatrix(V,totalColumn);
    MPI_Finalize();
}
