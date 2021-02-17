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
    //#pragma omp parallel for collapse(2)
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
    	//#pragma omp parallel for
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                for(int k = 0; k < size; k++)
                    r[i][j] = r[i][j] + a[i][k]*b[k][j];
}

double execution (float **a, float **b, float **r, int size, int threads){
    //omp_set_num_threads(8);
	//#pragma omp parallel for collapse(2)
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                r[i][j] = 0;
	//omp_set_num_threads(threads);
	double time;
    time=omp_get_wtime();
    multiply(a,b,r, size);
    time=omp_get_wtime()-time;
	//showMatrix(r, size);
    //cout << "\nExecution time: "<< time;
    return time;
}

void init(float **a, float **b, float **r, int size){
	//omp_set_num_threads(8);
	srand(time(NULL));
    int za, zb;
    double time=0;
    
    //#pragma omp parallel for
    for(int i = 0; i < size; i++){
        a[i] = (float *)malloc(size * sizeof(float));
        b[i] = (float *)malloc(size * sizeof(float));
        r[i] = (float *)malloc(size * sizeof(float));
    }
    
    double sec=omp_get_wtime();         
    /*#pragma omp sections
	{
		#pragma omp section 
		create_Matrix(a, size);
		
		#pragma omp section 
		create_Matrix(b, size);
	}
	*/
	create_Matrix(a, size);
	create_Matrix(b, size);
	sec= omp_get_wtime()-sec;
	cout<<"\ntempo setup: "<<sec;
	za=conta_zeri(a, size);
	zb=conta_zeri(b, size);
	int ntot=size*size;
	cout<<"\nNumero zeri di A: "<<za<<"\tNumero totale elementi di A: "<<ntot<<endl;
	cout<<"Numero zeri di B: "<<zb<<"\tNumero totale elementi di B: "<<ntot<<endl;
   /* cout << "\nMatrix A:\n";
    showMatrix(a, size);
    cout << "\nMatrix B:\n";
    showMatrix(b, size);
    cout << "\n\nA * B =\n";
   */
	
}

int main(){
	
	int dimension[] = { 3000, 3500, 4000};
	int threadcount[] = { 1};
	double avgtime, sum;
	ofstream outfile;
	outfile.open("Test_results_multiplication_seriale.txt");
	
	for (int i = 0; i < sizeof(dimension)/sizeof(dimension[0]); i++){
		float **a = (float **)malloc(dimension[i] * sizeof(float*));
    	float **b = (float **)malloc(dimension[i] * sizeof(float*));
    	float **r = (float **)malloc(dimension[i] * sizeof(float*));	
		cout <<"\n\n\nDimension: "<< dimension[i];
		outfile <<"\n\n\nDimension: "<< dimension[i];
		init(a, b, r, dimension[i]);
		for (int j = 0; j < sizeof(threadcount)/sizeof(threadcount[0]); j++){
			avgtime = 0;
			sum=0; 
			cout <<"\nNumber of threads: "<< threadcount[j];
			outfile <<"\nNumber of threads: "<< threadcount[j];
			for(int k=0;k<4;k++){	
				sum = sum+execution(a, b, r, dimension[i], threadcount[j]); 
			}
			avgtime =sum/4; 
			cout<<"\nTime: "<<avgtime<<"\n";
			outfile<<"\nTime: "<<avgtime<<"\n";
		}
		free(a);
    	free(b);
    	free(r);
	}
	outfile.close();
	return 0;
}


