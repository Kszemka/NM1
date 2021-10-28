////////////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt process
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define n 7
#define IT_MAX 12


void create_A(gsl_matrix * A);

gsl_vector* matrix_vector(gsl_vector * x, gsl_matrix * A);

gsl_vector* scalar_vector(double scalar, gsl_vector * x);

double do_scalar(gsl_vector * T, gsl_vector * x);

gsl_vector* vector_sub(gsl_vector * x, gsl_vector * v);

double do_norm(gsl_vector * x);

void vector_2_matrix(gsl_matrix * X, gsl_vector * x, int iter);

gsl_matrix* matrix_matrix(gsl_matrix * A, gsl_matrix * X);

gsl_matrix* matrix_transpos(gsl_matrix * A);

void matrix_fprintf(FILE * file, gsl_matrix * A);

gsl_vector* vector_from_matrix(gsl_matrix * X, int col);

////////////////////////////////////////////////////////////////////////////////

int main(){

	gsl_matrix* A = gsl_matrix_calloc(n, n);
	gsl_matrix* X = gsl_matrix_calloc(n, n);
	gsl_matrix* D = gsl_matrix_calloc(n, n);

	gsl_vector* x = gsl_vector_calloc(n);
	gsl_vector* lambda = gsl_vector_calloc(IT_MAX+1);
	gsl_vector* x_next = gsl_vector_calloc(n);

	create_A(A);

	FILE * lam = fopen("lambda.dat", "w");
	FILE * de = fopen("Macierz_D.dat", "w");

	gsl_vector_set_all(x, 1.0);

	for (int k = 0; k < n; k++){
		gsl_vector_set_all(x, 1.0);

		for (int i = 1; i <= IT_MAX; i++){

			gsl_vector_memcpy(x_next, matrix_vector(x, A));
			
			for (int j = 0; j < k; j++){
				
				gsl_vector_memcpy(x_next, vector_sub(x_next, scalar_vector(do_scalar(x_next, vector_from_matrix(X, j)), vector_from_matrix(X, j))));
			}

			gsl_vector_set(lambda, i, do_scalar(x_next, x)/do_scalar(x, x));
			fprintf(lam, "%d\t%15g\n", i, gsl_vector_get(lambda, i));
			gsl_vector_memcpy(x, scalar_vector(1./do_norm(x_next), x_next));
		}
		vector_2_matrix(X, x, k);
		fprintf(lam, "\n\n");
	}

	D = matrix_matrix(matrix_transpos(X), matrix_matrix(A, X));

	matrix_fprintf(de, D);

	gsl_vector_free(x);
	gsl_vector_free(x_next);
	gsl_vector_free(lambda);

	gsl_matrix_free(A);
	gsl_matrix_free(X);
	gsl_matrix_free(D);

	fclose(lam);
	fclose (de);
}


////////////////////////////////////////////////////////////////////////////////

void create_A(gsl_matrix * A){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			gsl_matrix_set(A, i, j, 1. / (sqrt(2.0 + abs(i - j))));
		}
	}
}

gsl_vector* matrix_vector(gsl_vector * x, gsl_matrix * A){
	double tmp = 0;
	gsl_vector* temp = gsl_vector_calloc(n);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			tmp += gsl_matrix_get(A, i, j) * gsl_vector_get(x, j);
		}

		gsl_vector_set(temp, i, tmp);
		tmp = 0;
	}

	return temp;
}

gsl_vector* scalar_vector(double scalar, gsl_vector * x){
	for (int i = 0; i < n; i++){
		gsl_vector_set(x, i, gsl_vector_get(x, i)*scalar);
	}

	return x;
}

double do_scalar(gsl_vector * T, gsl_vector * x){
	double tmp = 0;
	for (int i = 0; i < n; i++){
		tmp += gsl_vector_get(T, i)*gsl_vector_get(x, i);
	}
	return tmp;
}

gsl_vector* vector_sub(gsl_vector * x, gsl_vector * v){
	gsl_vector* tmp = gsl_vector_calloc(n);
	for (int i = 0; i < n; i++){
		gsl_vector_set(tmp, i, gsl_vector_get(x,i)-gsl_vector_get(v,i));
	}
	return tmp;
}

double do_norm(gsl_vector * x){
	return sqrt(do_scalar(x,x));
}


void vector_2_matrix(gsl_matrix * X, gsl_vector * x, int iter){
	for (int i = 0; i < n; i++){
		gsl_matrix_set(X, i, iter, gsl_vector_get(x, i));
	}
}

gsl_matrix* matrix_matrix(gsl_matrix * A, gsl_matrix * X){
	gsl_matrix* tmp = gsl_matrix_calloc(n, n);
	double temp;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			temp = 0;
			for (int k = 0; k < n; k++){
				temp += gsl_matrix_get(A,i,k) * gsl_matrix_get(X,k,j);
			}

			gsl_matrix_set(tmp, i, j, temp);
		}
	}
	return tmp;
}

gsl_matrix* matrix_transpos(gsl_matrix * A){
	gsl_matrix* tmp = gsl_matrix_calloc(n, n);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			gsl_matrix_set(tmp, i, j, gsl_matrix_get(A, j, i));
		}
	}
	return tmp;
}

void matrix_fprintf(FILE * file, gsl_matrix * A){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			fprintf(file, "%15g\t", gsl_matrix_get(A, i, j));
		}
		fprintf(file, "\n");
	}
}

gsl_vector* vector_from_matrix(gsl_matrix * X, int col){
	gsl_vector* tmp = gsl_vector_calloc(n);
	for (int i = 0; i < n; i++){
		gsl_vector_set(tmp, i, gsl_matrix_get(X, i, col));
	}
	return tmp;
}