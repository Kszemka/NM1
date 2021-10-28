////////////////////////////////////////////////////////////////////////////////////////
// Richardson extrapolation
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 8 
#define a 0
#define b 1

double function(double x);

void Richardson(double **arr);

void print(double **arr);

void f_print_array(double **array, FILE* file);

void Simpson(double **arr, FILE * file);

void Milne(double **arr, FILE * file);

///////////////////////////////////////////////////////////////////////////////////////////////

int main() {

    FILE* file1 = fopen("Simpson.dat", "w");
    FILE* file2 = fopen("Milne.dat", "w");

    double **arr = (double**) malloc((N+1)*sizeof(double*));
    for(int i=0; i<=N; i++){
        arr[i] = (double*) malloc((N+1)*sizeof(double));
    }

    Simpson(arr, file1);
    Milne(arr, file2);

    for(int i=0; i<N; i++){
        free(arr[i]);
    }
    free(arr);

    fclose(file1);
    fclose(file2);

}

////////////////////////////////////////////////////////////////////////////////////////////////
//// implementacja funkcji pomocniczych
///////////////////////////////////////////////////////////////////////////////////////////////

double function(double x){
    return (log(pow(x, 3) + 3*pow(x, 2) + x + 0.1)) * (sin(18*x));
}

void Richardson(double **arr){
    for (int i=1; i<=N; i++){
        for (int j=1; j<=i; j++){
            arr[i][j] = (pow(4, j)*arr[i][j-1] - arr[i-1][j-1]) / (pow(4, j) - 1);
        }
    }
}

void f_print_array(double **arr, FILE* file){  
    fprintf(file, "w\t D(w,0)\t\t D(w,w)\n");
    for (int i=0; i<=N; i++){
        fprintf(file, "%d\t %g\t %g\n", i, arr[i][0], arr[i][i]);
    }
}

void Simpson(double **arr, FILE * file){
    double suma;
    double hn;

    for(int i=0; i<=N; i++){
        suma = 0.0;
        hn = (b-a)/(pow(2,i+1));

        for(int j=0; j<=(pow(2,i+1)/2-1); j++){
            suma += (hn/3.0) * (function(a+2*j*hn) + 4*function(a+(2*j+1)*hn) + function(a+(2*j+2)*hn));
        }

        arr[i][0] = suma;
    }

    Richardson(arr);
    f_print_array(arr, file);
}    

void Milne(double **arr, FILE * file){
    double suma;
    double hn;

    for(int i=0; i<=N; i++){
        suma = 0.0;
        hn = (b-a)/pow(2,i+2);

        for(int j=0; j<=(pow(2,i+2)/4-1); j++){
            suma += ((4*hn)/90.0) * (7*function(a+4*j*hn) + 32*function(a+(4*j+1)*hn) + 12*function(a+(4*j+2)*hn) + 32*function(a+(4*j+3)*hn) + 7*function(a+(4*j+4)*hn));
        }
        arr[i][0] = suma;
    }

    Richardson(arr);
    f_print_array(arr, file);
}