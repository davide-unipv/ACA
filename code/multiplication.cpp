#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <omp.h>

using namespace std;

#define MAXNUMBER 100
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
        for(i = 0; i <size; i++)
            for(j = 0; j< size; j++)
                random[i][j] = rand() %(range);

}

int conta_zeri(float **matrix, int size){
	int n=0;
	for(int i=0; i <size;i++ ){
        for(int j=0; j< size; j++){
        	if(matrix[i][j]==0) n++;
		}
	}
	return n;
}


void multiply(float **a, float **b, float **r, int size){
    //the program will work only on square matrices
    //multiply a*b
    #pragma omp parallel for collapse(3)
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                for(int k = 0; k < size; k++)
                    r[i][j] = r[i][j] + a[i][k]*b[k][j];
}

double execution (int size, int threads){
    omp_set_num_threads(threads);
	srand(time(NULL));
    int i,j, za, zb;
    float **a = (float **)malloc(size * sizeof(float*));
    float **b = (float **)malloc(size * sizeof(float*));
    float **r = (float **)malloc(size * sizeof(float*));
    double time=0;
    
    #pragma omp parallel for
    for(int i = 0; i < size; i++){
        a[i] = (float *)malloc(size * sizeof(float));
        b[i] = (float *)malloc(size * sizeof(float));
        r[i] = (float *)malloc(size * sizeof(float));
    }
    #pragma omp parallel for collapse(2)
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                r[i][j] = 0;
                
    #pragma omp sections
	{
		#pragma omp section 
		create_Matrix(a, size);
		
		#pragma omp section 
		create_Matrix(b, size);
	}
	za=conta_zeri(a, size);
	zb=conta_zeri(b, size);
	int ntot=size*size;
	
	cout<<"\nNumero zeri di A: "<<za<<"\tNumero totale elementi di A: "<<ntot<<endl;
	cout<<"Numero zeri di B: "<<zb<<"\tNumero totale elementi di B: "<<ntot<<endl;
	
   /* cout << "\nMatrix A:\n";
    showMatrix(a, size);
    //create_Matrix(b, size);
    cout << "\nMatrix B:\n";
    showMatrix(b, size);
    cout << "\n\nA * B =\n";
    showMatrix(r, size);*/
    time=omp_get_wtime();
    multiply(a,b,r, size);
    time=omp_get_wtime()-time;
    //cout << "\nExecution time: "<< time;
    free(a);
    free(b);
    free(r);
    return time;
}

int main(){

	int dimension[] = { 500,1000,1500,2000,2500,3000 };
	int threadcount[] = { 1,2,4,6,8 };
	double avgtime;
	ofstream outfile;
	outfile.open("Test_results_multiplication.txt");

	for (int i = 0; i < sizeof(threadcount)/sizeof(threadcount[0]); i++)
	{
		cout <<"\n\nNumber of threads: "<< threadcount[i]<<"\n";
		outfile <<"\n\nNumber of threads: "<< threadcount[i]<<"\n";
		for (int j = 0; j < sizeof(dimension)/sizeof(dimension[0]); j++)
		{
			avgtime = 0; 
			cout <<"\nDimension: "<< dimension[j];
			outfile <<"\nDimension: "<< dimension[j];
			avgtime = execution(dimension[j], threadcount[i]); 
			cout<<"\nTime: "<<avgtime<<"\n";
			outfile<<"\nTime: "<<avgtime<<"\n";
		}
	}
	outfile.close();
	return 0;
}


