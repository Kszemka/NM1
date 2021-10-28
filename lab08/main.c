////////////////////////////////////////////////////////////////////////////////////////
// Polynomial interpolation using spline
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define min -5.0
#define max 5.0

void knot_normal(gsl_vector* x, int n, double step);

void fill_func1(gsl_vector* x, gsl_vector*  y, int n);

void fill_func2(gsl_vector* x, gsl_vector* y, int n);

void wyznacz_M(gsl_vector* x, gsl_vector* y, int n, gsl_vector* out, double alfa, double beta);

double wyznacz_Sx(gsl_vector *x, gsl_vector *y, gsl_vector *m, int n, double val);

double derivate(double x);

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {

	FILE* file1 = fopen("f1.dat", "w");
	FILE* file2 = fopen("f2.dat", "w");
	int array[] = {5,8,21}; 

	for(int list = 0; list <= 2 ; list++){

		int n = array[list];

		double step = 1.*(max-min)/(n-1);

		gsl_vector* x = gsl_vector_calloc(n+1); //tablica węzłów
		gsl_vector* m1 = gsl_vector_calloc(n+1); //tablica m funk1
		gsl_vector* m2 = gsl_vector_calloc(n+1); //tablica m funk2
		gsl_vector* y_1 = gsl_vector_calloc(n+1); //tablica wartosci funk1
		gsl_vector* y_2 = gsl_vector_calloc(n+1); //tablica wartosci funk2


		knot_normal(x, n, step);


		fill_func1(x, y_1, n); //wektor wartosci funkcji 1/(1+x^2)
		fill_func2(x, y_2, n); //wektor wartosci funkcji cos(2x)

		wyznacz_M(x, y_1, n, m1, 0, 0);
		wyznacz_M(x, y_2, n, m2, 0, 0);


		for(double k = min; k <= max; k += 0.01){
			if(n==10){
				fprintf(file1, "%.2f \t%g \t%g\n", k, wyznacz_Sx(x,y_1,m1,n,k), derivate(gsl_vector_get(x,k)));
				fprintf(file2, "%.2f \t%g\n", k, wyznacz_Sx(x,y_2,m2,n,k));		
			}
			else{
				fprintf(file1, "%.2f \t%g\n", k, wyznacz_Sx(x,y_1,m1,n,k));
				fprintf(file2, "%.2f \t%g\n", k, wyznacz_Sx(x,y_2,m2,n,k));	
			}	

		}

		gsl_vector_free(x);
		gsl_vector_free(y_1);
		gsl_vector_free(y_2);
		gsl_vector_free(m1);
		gsl_vector_free(m2);


		fprintf(file1, "\n\n");
		fprintf(file2, "\n\n");
	}	

	fclose(file1);
	fclose(file2);
}


////////////////////////////////////////////////////////////////////////////////////////////////

void knot_normal(gsl_vector * x, int n, double step){
	gsl_vector_set(x,0,min);
	gsl_vector_set(x,0,max);
	for(int i = 1; i <= n-1; i++){
		gsl_vector_set(x,i,min+step);
	}
}

void fill_func1(gsl_vector * x, gsl_vector * y, int n){
	for (int i = 0; i <= n; i++){
		gsl_vector_set(y,i,1/(1+pow(gsl_vector_get(x,i), 2)));
	}
}

void fill_func2(gsl_vector * x, gsl_vector * y, int n){
	for (int i = 0; i <= n; i++){
		gsl_vector_set(y,i,cos(2*gsl_vector_get(x,i)));
	}
}

void fill_Matrix(gsl_matrix* A, gsl_vector* lambda, gsl_vector* ni, int n){
	gsl_matrix_set(A, 0, 0, 1);
	gsl_matrix_set(A, n-1, n-1, 1);
	for(int i=2; i<n-1; i++){
		gsl_matrix_set(A, i, i, 2);
		gsl_matrix_set(A, i, i-1, gsl_vector_get(ni, i));
		gsl_matrix_set(A, i, i+1, gsl_vector_get(lambda, i));
	}
}

void wyznacz_M(gsl_vector* x, gsl_vector* y, int n, gsl_vector* m, double alfa, double beta){

	gsl_vector *d = gsl_vector_calloc(n+1);
	gsl_vector *lambda = gsl_vector_calloc(n);
	gsl_vector *ni = gsl_vector_calloc(n);
	gsl_matrix *A = gsl_matrix_calloc(n+1, n+1);

	gsl_vector_set(d, 0, alfa);
	gsl_vector_set(d, n - 1, beta);

	for (int i = 1; i < n - 1; i++){
		double h = gsl_vector_get(x, i) - gsl_vector_get(x,i-1);
		double h1 = gsl_vector_get(x, i+1) - gsl_vector_get(x,i);
		gsl_vector_set(lambda, i, h1/(h+h1));
		gsl_vector_set(ni, i, 1-gsl_vector_get(lambda, i));
		gsl_vector_set(d, i, (6.0 / (h+h1)) * ((gsl_vector_get(y,i + 1) - gsl_vector_get(y,i)) / h1) - (gsl_vector_get(y,i) - gsl_vector_get(y,i - 1) / h));
	}

	gsl_vector_set(d, 0, alfa);
	gsl_vector_set(d, n-1, beta);

	gsl_matrix_set_zero(A);
	fill_Matrix(A, lambda, ni, n);

	gsl_linalg_HH_svx(A, d);

	gsl_vector_memcpy(m, d);

	gsl_vector_free(d);
	gsl_vector_free(lambda);
	gsl_vector_free(ni);
	gsl_matrix_free(A);
}


double wyznacz_Sx(gsl_vector *x, gsl_vector *y, gsl_vector *m, int n, double val){
	double sx = 0;

	for (int i = 1; i < n; i++){
		double h = gsl_vector_get(x, i) - gsl_vector_get(x,i-1);
		if (val >= gsl_vector_get(x, i) && val <= gsl_vector_get(x, i+1)){
			sx += gsl_vector_get(m,i) * ((pow(gsl_vector_get(x, i + 1) - val, 3)) / (6 * h));
			sx += gsl_vector_get(m,i+1) * ((pow(val - gsl_vector_get(x,i), 3)) / (6 * h));
			sx += (val - gsl_vector_get(x,i)) * ((gsl_vector_get(y,i + 1) - gsl_vector_get(y, i) / h) - ((h / 6.0) * (gsl_vector_get(m,i + 1) - gsl_vector_get(m,i))));
			sx += gsl_vector_get(y,i) - gsl_vector_get(m,i) * (pow(h, 2) / 6.0);
		}
	}

	return sx;

}

double derivate(double x){
	double temp = 0.0;
	temp += 1/(1+pow((x - 0.01),2));
	temp -= 2 * 1/(1+pow(x,2));
	temp += 1/(1+pow((x + 0.01),2));
	temp /= pow(0.01, 2);

	return temp;
}