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
    srand(time(NULL));
        for(i = 0; i <size; i++)
            for(j = 0; j< size; j++)
                random[i][j] = rand() %(range);
}

void lu(float **a, float **l, float **u, int size){
    //int i = 0, j = 0, k = 0; //danno problemi queste variabili condivise??
    #pragma omp parallel for
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
    bool flag=false; //possibile problema di variabile condivisa
    #pragma omp parallel for // le iterazioni sono indipendenti tra di loro, quindi posso parallelizzare
    
	 for (int k = 0; k < size-1; k++){   
        int imax = 0;
        //foreach column i need to find which row has the maximum (in module) value
        for (int j = k; j < size; j++){
            //finding the maximum
            if (abs(a[j][k]) > abs(a[imax][k])){
                imax = j;
                //cout <<"\n iMax = " <<imax<<"riga:"<<j<<" colonna:"<<k;
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
  	//cout<< "\nForward colonna:"<<column<<"\n";
	//y[0] = p[0][column] / l[0][0];
	//cout<< "\n00 :"<<y[0]<<"\n";
	float sum=0;
    for(int i = 0; i< size; i++){
    	
        //y[i] = p[i][column];
        //cout <<"\ni: "<<i;
        for (int j = 0; j < i; j++){
            //the P[i] part has been done before
            sum = sum + l[i][j]*y[j] ;
        }
        y[i] = (p[i][column]-sum)/l[i][i];
        sum=0;
    }
}

void backwardSubst(float **u, float *y, float **a1, int column, int size){
    /**
     * solving A^-1 [i] = P[i]/sum of Lrow-i
    */
    //a1[size-1][column] = y[size-1]/u[size-1][size-1];
    //cout<< "\ncolumn "<<column<<"\nA1 "<<a1[size-1][column];
    //cout<< "\nColumn:"<<column;
    float sum;
    a1[size-1][column]=y[size-1]/ u[size-1][size-1];
    for(int i= size-2; i >= 0; i--){
        //a1[i][column] = y[i];
        //cout <<"\ni: "<<i;
		sum=y[i];
        for (int j = size-1; j > i; j--)
        {
            //cout<<"\tj: "<<j;
           sum = sum- u[i][j]*a1[j][column] ;
        }
        
        a1[i][column] = sum/u[i][i]; 
        sum=0;
    }
}

void findInverse(float **a, float **a1, float **l, float **u, float **p, int size){
    //i can do each column indipendentily
    #pragma omp parallel for 
    for(int i=0; i< size; i++){
        float* y = new float[size](); //vettore di puntatori tutto a 0. è di volta in volta la colonna che modifichiamo
        forwardSubst(l,p,i,y,size); //dal sito: p=b e y=d
        backwardSubst(u,y,a1,i,size);
        //a1 è l'inversa
    }

}

void multiply(float **a, float **b, float **r, int size){
    //the program will work only on square matrices
    //multiply a*b
    #pragma omp parallel for collapse(3) 
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                for(int k = 0; k < size; k++)
                    r[i][j] = r[i][j] + a[i][k]*b[k][j];
    //showMatrix(results);
}

double execution (int size,int threadcount){
    omp_set_num_threads(threadcount);

    int i, j;
	float **a = (float **)malloc(size * sizeof(float*));
    float **l = (float **)malloc(size * sizeof(float*));
    float **u = (float **)malloc(size * sizeof(float*));
    float **p = (float **)malloc(size * sizeof(float*));
    float **r = (float **)malloc(size * sizeof(float*));
    float **a1 = (float **)malloc(size * sizeof(float*));
    float **a_p = (float **)malloc(size * sizeof(float*));
    
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
    // create the identity pivot matrix
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
        
    create_Matrix(a,size);
    //cout << "\nMatrix A:\n";
    //showMatrix(a, size);
    //cout << "\nPivoting matrix:\n";
    //showMatrix(p,size);

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
    //cout << "\nPivoting....\n";
    pivoting(a_p,p,size);
    #pragma omp parallel for
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            u[i][j] = a_p[i][j];
            
    //cout << "\nMatrix A pivottata:\n";
    //showMatrix(a_p, size);
    //cout << "\nPivoting matrix:\n";
    //showMatrix(p, size);
    lu(a_p,l,u,size); //come u posso passargli a_p
    /* TEST LU
    cout << "\nL matrix:\n";
    showMatrix(l, size);
    cout << "\nU matrix:\n";
    showMatrix(u, size);
    cout << "\nL*U=R";
    multiply(l,u,r, size);
    showMatrix(r, size);
    */
    
   /**
    * solving the system LUxi=Pi where i is the column and x is the column
    * of the inverse to find
    */
    //cout << "\nFinding the inverse...";
    findInverse(a_p,a1,l,u,p,size);
    
    //cout<< "\nMatrix a^-1:\n";
    //showMatrix(a1, size);

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
    int dimension[] = { 500,1000,1500,2000,2500,3000 };
	int threadcount[] = { 5,6,8 };
    double avgtime;
	ofstream outfile;
	outfile.open("Test_results_inverse.txt");
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
