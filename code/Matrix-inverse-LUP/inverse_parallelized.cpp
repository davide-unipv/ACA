#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <conio.h>
#include <omp.h>

using namespace std;
#define SIZE 3000
#define THREADNUMB 24
#define MAXNUMBER 50
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
    //srand(static_cast<unsigned>(time(0))+number); 
    srand(time(NULL));
    //srand(0);
    //#pragma omp parallel for collapse(2)
        for(i = 0; i <SIZE; i++)
            for(j = 0; j< SIZE; j++)
                random[i][j] = rand() %(range)+i*j;
}

void lu(float **a, float **l, float **u){
    int i = 0, j = 0, k = 0;
    
    for (k = 0; k < SIZE; k++){
        for (i = k+1; i < SIZE; i++){
                l[i][k] = u[i][k] / u[k][k];
                for (j = k; j < SIZE; j++){
                    u[i][j] = u[i][j] - l[i][k] * u[k][j];
                }
            }
        }
}

void pivoting(float **a, float **p){
    
    //#pragma omp parallel for //performance decreases
    for (int k = 0; k < SIZE-1; k++){   
        int imax = 0;
        //foreach column i need to find which row has the maximum (in module) value
        //#pragma omp parallel for
            for (int j = k; j < SIZE; j++){
            //finding the maximum
            if (abs(a[j][k]) > abs(a[imax][k])){
                imax = j;
                //cout <<"\nNew iMax = " << imax;
            }
        }
        //now i swap the row imax and k
        //cout <<"\nSwapping row "<<k << "and "<< imax;

        swap(a[k],a[imax]);
        swap(p[k],p[imax]);
    }
}

void forwardSubst(float **l, float **p, float *y, int column){
    /**
     * solving yi=Pi/ sum of Lrow-i
    */
    //y[0]=p[0][column] /l[0][0];
    //cout<< "\nColumn:"<<column;

	y[0] = p[0][column] / l[0][0];
    for(int i = 1; i< SIZE; i++){
        y[i] = p[i][column];
        //cout <<"\ni: "<<i;
        for (int j = 0; j < i; j++)
        {
            //cout<<"\tj: "<<j;
            //the P[i] part has been done before
            y[i] = y[i] - l[i][j]*y[j] ;
        }
        y[i] = y[i]/l[i][i];
    }
}

void backwardSubst(float **u, float *y, float **a1, int column){
    /**
     * solving A^-1 [i] = P[i]/sum of Lrow-i
    */
    //a1[SIZE-1][column] = y[SIZE-1]/u[SIZE-1][SIZE-1];
    //cout<< "\ncolumn "<<column<<"\nA1 "<<a1[SIZE-1][column];
    //cout<< "\nColumn:"<<column;
    a1[SIZE-1][column] = y[SIZE-1] / u[SIZE-1][SIZE-1];
    for(int i= SIZE-2; i >= 0; i--){
        //a1[i][column] = y[i];
        //cout <<"\ni: "<<i;

        for (int j = SIZE-1; j > i; j--)
        {
            //cout<<"\tj: "<<j;
            a1[i][column] = a1[i][column] - u[i][j]*a1[j][column] ;
        }
        a1[i][column] = a1[i][column]/u[i][i]; 
    }
}

void findInverse(float **a, float **a1, float **l, float **u, float **p){
    /**
     * foreach column i solve the system LUai=Pi with the i-th column
     * by using the forward substitution method
     * LYi=Pi
     * Uxi=Yi
     */

    //i can only parallelize this, so i can do each column indipendentily
    #pragma omp parallel for 
    for(int i=0; i< SIZE; i++){
        //float *y = (float *)malloc(SIZE * sizeof(float*));
        //y = yi;
        float* y = new float[SIZE]();
        forwardSubst(l,p,y,i);
        backwardSubst(u,y,a1,i);
    }

}

void multiply(float **a, float **b, float **r){
    //the program will work only on square matrices
    //multiply a*b
    #pragma omp parallel for 
        for(int i = 0; i < SIZE; i++)
            for(int j = 0; j < SIZE; j++)
                for(int k = 0; k < SIZE; k++)
                    r[i][j] = r[i][j] + a[i][k]*b[k][j];
    //showMatrix(results);
}


int main(void){
    omp_set_num_threads(THREADNUMB);

    int i, j;
	float **a = (float **)malloc(SIZE * sizeof(float*));
    float **l = (float **)malloc(SIZE * sizeof(float*));
    float **u = (float **)malloc(SIZE * sizeof(float*));
    float **p = (float **)malloc(SIZE * sizeof(float*));
    float **r = (float **)malloc(SIZE * sizeof(float*));
    float **a1 = (float **)malloc(SIZE * sizeof(float*));
    //float *y = (float *)malloc(SIZE * sizeof(float*)); 
    
    float **a_p = (float **)malloc(SIZE * sizeof(float*));
    
    #pragma omp parallel for 
    for(int i = 0; i < SIZE; i++){
        a[i] = (float *)malloc(SIZE * sizeof(float));
        a1[i] = (float *)malloc(SIZE * sizeof(float));
        l[i] = (float *)malloc(SIZE * sizeof(float));
        u[i] = (float *)malloc(SIZE * sizeof(float));
        p[i] = (float *)malloc(SIZE * sizeof(float));
        r[i] = (float *)malloc(SIZE * sizeof(float));
        a_p[i] = (float *)malloc(SIZE * sizeof(float));
    }
    // create the identity pivot matrix
    #pragma omp parallel for
        for(int i = 0; i < SIZE; i++) {
            p[i][i] = 1;
            for(int j = 0; j < SIZE; j++) {
                a1[i][j] = 0;
                r[i][j] = 0;

                if (i != j) {
                    p[i][j] = 0;
                }
            }
        }

    
    
    create_Matrix(a,10);
    cout << "\nMatrix A:\n";
    //showMatrix(a);
    //cout << "\nPivoting matrix:\n";
    //showMatrix(p);
    double time= omp_get_wtime();
    
    #pragma omp parallel for   
    for (i = 0; i < SIZE; i++) {
        l[i][i] = 1;
        for (j = 0; j < SIZE; j++) {
            a_p[i][j] = a[i][j];

            if (i != j)
                l[i][j] = 0;
        }
    }
    
    /**
     * P*a = l*u and LUx =p
     * where p is the column to pivot of the b matrix and x is a column of the inverse
    */
    //cout << "\nPivoting....\n";
    pivoting(a_p,p);
    
    #pragma omp parallel for
    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++)
            u[i][j] = a_p[i][j];

    /**cout << "\nMatrix A:\n";
    showMatrix(a);
    cout << "\nPivoting matrix:\n";
    showMatrix(p);
    */
    lu(a_p,l,u);
     
    /* TEST LU
    // /**
    cout << "\nL matrix:\n";
    showMatrix(l);
    cout << "\nU matrix:\n";
    showMatrix(u);

    cout << "\nL*U";
    multiply(l,u,r);
    showMatrix(r);
    
    
    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++)
            r[i][j] = 0;
    */

    // */

   /**
    * solving the system LUxi=Pi where i is the column and x is the column
    * of the inverse to find
    */
    cout << "\nFinding the inverse...";
    findInverse(a_p,a1,l,u,p);
    
    cout<< "\nMatrix a^-1:\n";
    //showMatrix(a1);

    time = omp_get_wtime()-time;
    cout << "\nExecution time: "<< time << "\n\n";
    
    //multiply(a,a1,r);
    //showMatrix(r);
    getch();
    return 0;
}