////////////////////////////////////////////////////////////////////////////////////////
// Polynomial regression
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define n 11 //101
#define sigma 4
#define x0 2
#define max (3.*sigma + (double)x0)
#define min (-3.*sigma + (double)x0)
#define a0 (-1.*pow(x0,2)/2./pow(sigma,2))
#define a1 ((double)x0/pow(sigma,2))
#define a2 (-1./2./pow(sigma,2))

void set_vector(gsl_vector* v, const double step); 

void set_matrix(gsl_matrix* M, const double step);

void vector_printf(gsl_vector * v);

void matrix_printf(gsl_matrix * M);

double aprox_function(int value, const double step);

double x_value(int value, const double step);

double U ();

////////////////////////////////////////////////////////////////////////////////////////

int main() {

	srand(time(NULL));

	FILE* file1 = fopen("G.dat", "w");
	FILE* file2 = fopen("pkt.dat", "w");

	double step = 1.*(max-min)/(n-1);

	gsl_vector* r = gsl_vector_calloc(sigma); 
	gsl_matrix* G = gsl_matrix_calloc(sigma, sigma); 

	set_vector(r, step);
	set_matrix(G, step);

	printf("Wektor r:\n");
	vector_printf(r);
	printf("Macierz G:\n");
	matrix_printf(G);

	gsl_linalg_HH_svx(G,r);

	for (int i = 1; i <= n; i++){
        fprintf(file1," %f \t%f\n", x_value(i, step), exp(gsl_vector_get(r, 0) + gsl_vector_get(r, 1) * x_value(i, step) + gsl_vector_get(r, 2) * pow(x_value(i, step), 2) + gsl_vector_get(r, 3) * pow(x_value(i, step), 3)));
        fprintf(file2," %f \t%f\n", x_value(i, step), exp(aprox_function(i, step)));
    }

    fprintf(file1, "\n\n");
    fprintf(file2, "\n\n");

    for (int i = 1; i <= n; i++){
        fprintf(file1," %f \t%f\n", x_value(i, step), exp(gsl_vector_get(r, 0) + gsl_vector_get(r, 1) * x_value(i, step) + gsl_vector_get(r, 2) * pow(x_value(i, step), 2) + gsl_vector_get(r, 3) * pow(x_value(i, step), 3)));
        fprintf(file2," %f \t%f\n", x_value(i, step), exp(aprox_function(i, step))*(1.+0.5*(U()-0.5)));
    }

    fprintf(file1, "\n");
    fprintf(file2, "\n");

	gsl_vector_free(r);
	gsl_matrix_free(G);

	fclose(file1);
	fclose(file2);
}


////////////////////////////////////////////////////////////////////////////////////////////////

void set_vector(gsl_vector* v, const double step){
    for (int i = 1; i <= sigma; i++){
        double value = 0.;
        for (int j = 1; j <= n; j++){
            value += aprox_function(j, step) * pow(x_value(j, step), i - 1);
        }
        gsl_vector_set(v, i - 1, value);
    }
}

void set_matrix(gsl_matrix* M, const double step){
    for (int i = 1; i <= sigma; i++){
        for (int j = 1; j <= sigma; j++){
            double value = 0.;
            for (int k = 1; k <= n; k++){
                value += pow(x_value(k, step), i + j - 2);
            }
            gsl_matrix_set(M, i - 1, j - 1, value);
        }
    }
}

void vector_printf(gsl_vector * v){
	for (int i = 0; i < sigma; i++){
		printf("%f\t", gsl_vector_get(v, i));
	}
	printf("\n");
}

void matrix_printf(gsl_matrix * M){
	for (int i = 0; i < sigma; i++){
		for (int j = 0; j < sigma; j++){
			printf("%f\t", gsl_matrix_get(M, i, j));
		}
		printf("\n");
	}
}

double aprox_function(int value, const double step){    
    int temp = 0;
    for (double x = min; x <= max + step; x += step){
        temp++;
        if (temp == value){
            return a0 + (a1 * x) + (a2 * x * x);
        }
    }
    return 0;
}

double x_value(int value, const double step){
    int temp = 0;
    for (double x = min; x <= max + step; x += step){
        temp++;
        if (temp == value){
            return x;
        }
    }
    return 0;
}

double U (){
 return rand()/(RAND_MAX+1.0);
}
