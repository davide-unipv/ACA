#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <omp.h>

using namespace std;
#define MAXNUMBER 50
#define MINNUMBER 0
 
void showMatrix(float **matrix, int size){
    cout << "\n";
    for(int i=0; i <size;i++ ){
        for(int j=0; j< size; j++){
            cout << matrix[i][j];
            cout << "\t";
        }
        cout << "\n";
    }
}

void create_Matrix (float **random, int size){
    
    int i, j;
    int range = MAXNUMBER - MINNUMBER;
    //srand(static_cast<unsigned>(time(0))+number); 
    srand(time(NULL));
    //srand(0);
        for(i = 0; i <size; i++)
            for(j = 0; j< size; j++)
                random[i][j] = rand() %(range);
}

void lu(float **a, float **l, float **u, int size){
    int i = 0, j = 0, k = 0;
    
    for (k = 0; k < size; k++){
        for (i = k+1; i < size; i++){
                l[i][k] = u[i][k] / u[k][k];
                for (j = k; j < size; j++){
                    u[i][j] = u[i][j] - l[i][k] * u[k][j];
                }
            }
        }
}

void pivoting(float **a, float **p, int size){
    
    //#pragma omp parallel for //performance decreases because of the memory conflict that if not managed produces a wrong output
    for (int k = 0; k < size-1; k++){   
        int imax = 0;
        //foreach column i need to find which row has the maximum (in module) value
            for (int j = k; j < size; j++){
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

void forwardSubst(float **l, float **p, float *y, int column, int size){
    /**
     * solving yi=Pi/ sum of Lrow-i
    */
    //y[0]=p[0][column] /l[0][0];
    //cout<< "\nColumn:"<<column;

	y[0] = p[0][column] / l[0][0];
    for(int i = 1; i< size; i++){
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

void backwardSubst(float **u, float *y, float **a1, int column, int size){
    /**
     * solving A^-1 [i] = P[i]/sum of Lrow-i
    */
    //a1[size-1][column] = y[size-1]/u[size-1][size-1];
    //cout<< "\ncolumn "<<column<<"\nA1 "<<a1[size-1][column];
    //cout<< "\nColumn:"<<column;
    a1[size-1][column] = y[size-1] / u[size-1][size-1];
    for(int i= size-2; i >= 0; i--){
        //a1[i][column] = y[i];
        //cout <<"\ni: "<<i;

        for (int j = size-1; j > i; j--)
        {
            //cout<<"\tj: "<<j;
            a1[i][column] = a1[i][column] - u[i][j]*a1[j][column] ;
        }
        a1[i][column] = a1[i][column]/u[i][i]; 
    }
}

void findInverse(float **a, float **a1, float **l, float **u, float **p, int size){
    /**
     * foreach column i solve the system LUai=Pi with the i-th column
     * by using the forward substitution method
     * LYi=Pi
     * Uxi=Yi
     */

    //i can only parallelize this, so i can do each column indipendentily
    #pragma omp parallel for 
    for(int i=0; i< size; i++){
        float* y = new float[size]();
        forwardSubst(l,p,y,i,size);
        backwardSubst(u,y,a1,i,size);
    }

}

double execution (int size,int threadcount){
    omp_set_num_threads(threadcount);

    int i, j;
	float **a = (float **)malloc(size * sizeof(float*));
    float **l = (float **)malloc(size * sizeof(float*));
    float **u = (float **)malloc(size * sizeof(float*));
    float **p = (float **)malloc(size * sizeof(float*));
    //float **r = (float **)malloc(size * sizeof(float*));
    float **a1 = (float **)malloc(size * sizeof(float*));
    //float *y = (float *)malloc(size * sizeof(float*)); 
    
    float **a_p = (float **)malloc(size * sizeof(float*));
    
    #pragma omp parallel for 
    for(int i = 0; i < size; i++){
        a[i] = (float *)malloc(size * sizeof(float));
        a1[i] = (float *)malloc(size * sizeof(float));
        l[i] = (float *)malloc(size * sizeof(float));
        u[i] = (float *)malloc(size * sizeof(float));
        p[i] = (float *)malloc(size * sizeof(float));
        //r[i] = (float *)malloc(size * sizeof(float));
        a_p[i] = (float *)malloc(size * sizeof(float));
    }
    // create the identity pivot matrix
    #pragma omp parallel for
        for(int i = 0; i < size; i++) {
            p[i][i] = 1;
            for(int j = 0; j < size; j++) {
                a1[i][j] = 0;
                //r[i][j] = 0;

                if (i != j) {
                    p[i][j] = 0;
                }
            }
        }

    
    
    create_Matrix(a,size);
    //cout << "\nMatrix A:\n";
    //showMatrix(a);
    //cout << "\nPivoting matrix:\n";
    //showMatrix(p);
    double time= omp_get_wtime();
    
    #pragma omp parallel for   
    for (i = 0; i < size; i++) {
        l[i][i] = 1;
        for (j = 0; j < size; j++) {
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
    pivoting(a_p,p,size);
    
    #pragma omp parallel for
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            u[i][j] = a_p[i][j];

    /**cout << "\nMatrix A:\n";
    showMatrix(a);
    cout << "\nPivoting matrix:\n";
    showMatrix(p);
    */
    lu(a_p,l,u,size);
     
    /* TEST LU
    // /**
    cout << "\nL matrix:\n";
    showMatrix(l);
    cout << "\nU matrix:\n";
    showMatrix(u);

    cout << "\nL*U";
    multiply(l,u,r);
    showMatrix(r);
    
    
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            r[i][j] = 0;
    */

    // */

   /**
    * solving the system LUxi=Pi where i is the column and x is the column
    * of the inverse to find
    */
    //cout << "\nFinding the inverse...";
    findInverse(a_p,a1,l,u,p,size);
    
    //cout<< "\nMatrix a^-1:\n";
    //showMatrix(a1);

    time = omp_get_wtime()-time;
    //cout << "\nExecution time: "<< time << "\n\n";
    
    //multiply(a,a1,r);
    //showMatrix(r);
    free(*a);
    free(*a_p);
    free(*a1);
    free(*l);
    free(*u);
    free(*p);
    return time;
}

int main(int argc,char **argv){
        

    int dimension[]={2000,4000,6000,8000,10000,20000};
    int threadcount[]={1,2,4,6,8,12,24};
    double avgtime;
	ofstream outfile;
	outfile.open("Test_results_inverse.txt");
    //int size = atoi(argv[1]);
    for (int i = 0; i < 7; i++)
    {
        cout<<"\n\nThreads: "<<threadcount[i]<<"\nSize:\tTime AVG:\n";
		outfile <<"\n\nThreads: " << threadcount[i] << "\nSize:\tTime AVG:\n";
        for (int j = 0; j < 6; j++)
        {
            avgtime = 0; //reinitilize it
            cout<<dimension[j]<<"\t";
			outfile << dimension[j] << "\t";
            //for (int k =1; k <= 5; k++){
				avgtime = execution(dimension[j], threadcount[i]);//avgtime + execution(dimension[j],threadcount[i]);
            //}
            //avgtime = avgtime/5.0F;
            cout<<avgtime<<"\n";
			outfile << avgtime << "\n";
        }
    }
	outfile.close();
    return 0;
}