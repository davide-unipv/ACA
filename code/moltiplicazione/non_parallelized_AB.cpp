#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <conio.h>
#include <omp.h>

using namespace std;
#define SIZE 2000
#define MAXNUMBER 10
#define MINNUMBER 0
 
void showMatrix(float **matrix){
    cout << "\n";
    for(int i=0; i <SIZE;i++ ){
        for(int j=0; j< SIZE; j++){
            cout << matrix[i][j];
            cout << "\t";
        }
        cout << "\n";
    }
}

void create_Matrix (float **random,int number){
    
    int i, j;
    int range = MAXNUMBER - MINNUMBER;
    srand(static_cast<unsigned>(time(0))+number); 
    //srand(time(NULL));
    #pragma omp parallel for
        for(i = 0; i <SIZE; i++)
            for(j = 0; j< SIZE; j++)
                random[i][j] = rand() %(range);
}

void multiply(float **a, float **b, float **r){
    //the program will work only on square matrices
    //multiply a*b
    for(int i = 0; i < SIZE; i++)
        for(int j = 0; j < SIZE; j++)
            for(int k = 0; k < SIZE; k++)
                r[i][j] = r[i][j] + a[i][k]*b[k][j];
    //showMatrix(results);
}

int main(void){
    int i,j,n = SIZE;
	//float a[SIZE][SIZE], b[SIZE][SIZE], at[SIZE][SIZE], bt[SIZE][SIZE];
    float **a = (float **)malloc(SIZE * sizeof(float*));
    float **b = (float **)malloc(SIZE * sizeof(float*));
    float **r = (float **)malloc(SIZE * sizeof(float*));
    
    for(int i = 0; i < SIZE; i++){
        a[i] = (float *)malloc(SIZE * sizeof(float));
        b[i] = (float *)malloc(SIZE * sizeof(float));
        r[i] = (float *)malloc(SIZE * sizeof(float));
    }
    #pragma omp parallel for 
        for(int i = 0; i < SIZE; i++)
            for(int j = 0; j < SIZE; j++)
                r[i][j] = 0;
    
    
    create_Matrix(a,1);
    //cout << "\nMatrix A:\n";
    //showMatrix(a);
    create_Matrix(b,2);
    //cout << "\nMatrix B:\n";
    //showMatrix(b);
    cout << "\n\nA * B =\n";
    double time= omp_get_wtime();
    multiply(a,b,r);
    time = omp_get_wtime()-time;
    cout << "\nExecution time: "<< time;
    getch();
    return 0;
}