////////////////////////////////////////////////////////////////////////////////////////
// Polynomial interpolation
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min -5
#define max 5

double interpolation(double* x, double* y, double value, int n);

void fill_func(double* x, double* y, int n);

void knot_normal(double* x, int n, double step);

void knot_Cz(double* x, int n);

//////////////////////////////////////////////////////////////////////////////////

int main() {

	FILE* file1 = fopen("zad_1.dat", "w");
	FILE* file2 = fopen("zad_2.dat", "w");

	for(int n = 5; n <= 20; n += 5){

		double step = 1.*(max-min)/n;

		double* x_n = malloc((n+1) * sizeof(double)); //tablica węzłów
		double* x_C = malloc((n+1) * sizeof(double)); //tablica węzłów Czebyszew
		double* y_n = malloc((n+1) * sizeof(double)); //tablica wartosci dla x_n
		double* y_C = malloc((n+1) * sizeof(double)); //tablica wartosci dla x_C


		knot_normal(x_n, n, step);
		knot_Cz(x_C, n);

		fill_func(x_n, y_n, n); //wektor wartosci funkcji dla wezlow rozklad normalny
		fill_func(x_C, y_C, n); //wektor wartosci funkcji dla wezlow Czebyszew

		for(double k = -5.0; k <= 5.0; k += 0.01){
			fprintf(file1, "%.2f \t%g\n", k, interpolation(x_n,y_n,k,n));
			fprintf(file2, "%.2f \t%g\n", k, interpolation(x_C,y_C,k,n));
		}

		free(x_n);
		free(x_C);
		free(y_n);
		free(y_C);

		fprintf(file1, "\n\n");
		fprintf(file2, "\n\n");
		printf("\n\n");

	}	

	fclose(file1);
	fclose(file2);
}


///////////////////////////////////////////////////////////////////////////////////////////


double interpolation(double* x, double* y, double value, int n) {
	double result = 0.0;
	double licznik;
	double mianownik;

	for (int i = 0; i <= n; i++){
		licznik = 1.0;
		mianownik = 1.0;
		for (int j = 0; j <=n ; j++){
			if (i != j){
				mianownik *= (x[i] - x[j]);
				licznik *= (value - x[j]);
			}	
		}

		result += y[i] * (licznik / mianownik);
	}

	return result;
}

void fill_func(double *x, double *y, int n){
	for (int i = 0; i <= n; i++){
		y[i] = exp (-pow(x[i], 2));
	}
}

void knot_normal(double * x, int n, double step){
	x[0] = min;
	x[n] = max;
	for(int i = 1; i <= n-1; i++){
		x[i] = min + (i*step);
	}
}

void knot_Cz(double* x, int n){
	for (int i = 0; i <= n; i++){
		x[i] = 0.5*((max-min)*(cos(M_PI*((double)((2*i)+1)/((2*n)+2)))+(min+max)));
	}
}