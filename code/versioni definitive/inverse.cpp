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
    srand(time(NULL));
    for(i = 0; i <size; i++)
        for(j = 0; j< size; j++)
            random[i][j] = rand() %(range);
}

void lu(float **a, float **l, float **u, int size){
    for (int k = 0; k < size; k++){
        for (int i = k+1; i < size; i++){
            l[i][k] = u[i][k] / u[k][k];
            for (int j = k; j < size; j++){
                u[i][j] = u[i][j] - l[i][k] * u[k][j];
        	}
        }
    }
}

void pivoting(float **a, float **p, int size){
    bool flag=false; 
	for (int k = 0; k < size-1; k++){   //k=colonna
    	int imax = k;
        for (int j = k; j < size; j++){ //j=riga
            if (a[j][k] > a[imax][k]){ //trovo l'indice massimo
				imax = j;
                flag=true; 	
            }
        }
        if(flag){
			//cout <<"\nSwapping row "<<k << "and "<< imax;
        	swap(a[k],a[imax]);
        	swap(p[k],p[imax]);
        	flag=false;
    	}
	}
}

void forwardSubst(float **l, float **p,int column, float *y, int size){
	double sum=0;
    for(int i = 0; i< size; i++){
        for (int j = 0; j < i; j++){
            sum = sum + l[i][j]*y[j] ;
        }
        y[i] = (p[i][column]-sum)/l[i][i];
        sum=0;
    }
}

void backwardSubst(float **u, float *y, float **a1, int column, int size){
    double sum;
    a1[size-1][column]=y[size-1]/ u[size-1][size-1];
    for(int i= size-2; i >= 0; i--){
		sum=y[i];
        for (int j = size-1; j > i; j--){
           sum = sum- u[i][j]*a1[j][column] ;
        }
        a1[i][column] = sum/u[i][i]; 
        sum=0;
    }
}

void findInverse(float **a, float **a1, float **l, float **u, float **p, int size){
    #pragma omp parallel for  //i can do each column indipendentily
    for(int i=0; i< size; i++){
        float* y = new float[size](); 
        forwardSubst(l,p,i,y,size); 
        backwardSubst(u,y,a1,i,size);
    }
}

void multiply(float **a, float **b, float **r, int size){ //per testare l'inversa
    #pragma omp parallel for 
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                for(int k = 0; k < size; k++)
                    r[i][j] = r[i][j] + a[i][k]*b[k][j];
}

int conta_zeri(float **matrix, int size){ //conta il numero di zeri nella matrice
	int n=0;
	for(int i=0; i <size;i++ ){
        for(int j=0; j< size; j++){
        	if(matrix[i][j]==0) n++;
		}
	}
	return n;
}

double execution (float **a, float **l, float **u, float **p, float **r, float **a1, float **a_p, int size, int threadcount){
	omp_set_num_threads(8);
	#pragma omp parallel for
    for(int i = 0; i < size; i++) {
        p[i][i] = 1;
        for(int j = 0; j < size; j++) {
            a1[i][j] = 0;
            r[i][j] = 0;
            if (i != j) {
                p[i][j] = 0;
            }
        }
    }
    omp_set_num_threads(threadcount);
    
	double time_exec= omp_get_wtime();
	#pragma omp parallel for   
    for (int i = 0; i < size; i++) {
        l[i][i] = 1;
        for (int j = 0; j < size; j++) {
            a_p[i][j] = a[i][j];
            if (i != j)
                l[i][j] = 0;
        }
    }
    //double pivot=omp_get_wtime();
    pivoting(a_p,p,size);
	/*pivot=omp_get_wtime()-pivot;
	cout<<"\ntempo pivoting: "<<pivot; */
	    
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            u[i][j] = a_p[i][j];
    
    /*cout << "\nMatrix A pivottata:\n";
    showMatrix(a_p, size);*/
    
    double lu_time=omp_get_wtime();
    lu(a_p,l,u,size); 
	lu_time=omp_get_wtime()-lu_time;
    cout<<"\ntempo lu: "<<lu_time;        
    /*TEST LU 
    cout << "\nL matrix:\n";
    showMatrix(l, size);
    cout << "\nU matrix:\n";
    showMatrix(u, size); 
    cout << "\nL*U=R";
    multiply(l,u,r, size);
    showMatrix(r, size);
	cout << "\nFinding the inverse...";
   */
    findInverse(a_p,a1,l,u,p,size);
    //cout<< "\nMatrix a^-1:\n";
    //showMatrix(a1, size);
    time_exec = omp_get_wtime()-time_exec;
    /*cout << "\nTempo inversa: "<< inv;
    cout << "\nExecution time: "<< time << "\n\n";
    cout<<"\nmoltiplicazione a*a-1: "<<endl;
    multiply(a,a1,r, size);
    showMatrix(r, size);*/
    return time_exec;
}
void init(float **a, float **l, float **u, float **p, float **r, float **a1, float **a_p, int size){
	double time_setup= omp_get_wtime();
	omp_set_num_threads(8);
	#pragma omp parallel for 
    for(int i = 0; i < size; i++){
        a[i] = (float *)malloc(size * sizeof(float));
        a1[i] = (float *)malloc(size * sizeof(float));
        l[i] = (float *)malloc(size * sizeof(float));
        u[i] = (float *)malloc(size * sizeof(float));
        p[i] = (float *)malloc(size * sizeof(float));
        r[i] = (float *)malloc(size * sizeof(float));
        a_p[i] = (float *)malloc(size * sizeof(float));
    }
	
    create_Matrix(a,size);	
	
	time_setup = omp_get_wtime()-time_setup;
	cout<<"\ntempo setup: "<<time_setup;
	int za=conta_zeri(a, size);
	int ntot=size*size;
	cout<<"\nNumero zeri di A: "<<za<<"\tNumero totale elementi di A: "<<ntot<<endl;

    /*cout << "\nMatrix A:\n";
    showMatrix(a, size);*/
}

void free_mem(float **a, float **l, float **u, float **p, float **r, float **a1, float **a_p){
	free(*a);
    free(*a_p);
    free(*a1);
    free(*l);
    free(*u);
    free(*p);
    free(*r);
}

int main(int argc,char **argv){
    int dimension[] = {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000};
	int threadcount[] = {2, 4, 8, 12, 16, 20, 24};
    double avgtime, sum;
	ofstream outfile;
	outfile.open("Test_results_inverse.txt");
	for (int i = 0; i < sizeof(dimension)/sizeof(dimension[0]); i++){
		float **a = (float **)malloc(dimension[i] * sizeof(float*));
    	float **l = (float **)malloc(dimension[i] * sizeof(float*));
    	float **u = (float **)malloc(dimension[i] * sizeof(float*));
   		float **p = (float **)malloc(dimension[i] * sizeof(float*));
    	float **r = (float **)malloc(dimension[i] * sizeof(float*));
    	float **a1 = (float **)malloc(dimension[i] * sizeof(float*));
    	float **a_p = (float **)malloc(dimension[i] * sizeof(float*));
		cout <<"\n\n\nDimension: "<< dimension[i];
		outfile <<"\n\n\nDimension: "<< dimension[i];
		init(a, l, u, p, r, a1, a_p, dimension[i]);
		for (int j = 0; j < sizeof(threadcount)/sizeof(threadcount[0]); j++){
			avgtime = 0;
			sum=0; 
			cout <<"\nNumber of threads: "<< threadcount[j];
			outfile <<"\nNumber of threads: "<< threadcount[j];
			for(int k=0; k<4; k++){
				sum = sum+execution(a, l, u, p, r, a1, a_p, dimension[i], threadcount[j]); 	
			}
			avgtime=sum/4;
			cout<<"\nTime: "<<avgtime<<"\n";
			outfile<<"\nTime: "<<avgtime<<"\n";
		}
		free_mem(a, l, u, p, r, a1, a_p);
	}
	outfile.close();
    return 0;
}
