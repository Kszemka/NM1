////////////////////////////////////////////////////////////////////////////////////////
// Newton's method 
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5
#define IT_MAX 30

double licz_r(double * a, double * b, int n, double xj) {
	b[n] = 0;
	for(int i = n-1; i >= 0; i--) b[i] = a[i+1] + xj*b[i+1];
	return (a[0] + xj*b[0]);
}

///////////////////////////////////////////////////////////////////////////////////

int main() {

	double* a = malloc((N+1)*sizeof(double));
	double* b = malloc((N+1)*sizeof(double));
	double* c = malloc((N+1)*sizeof(double));

	FILE* file = fopen("out.dat", "w");

	fprintf(file, "%s \t| %s \t| %s \t| %s \t| %s \t|\n\n", "L", "it", "x0", "R_j", "R_j'");

	a[0] = 240.;
	a[1] = -196.;
	a[2] = -92.;
	a[3] = 33.;
	a[4] = 14.;
	a[5] = 1.;

	double x0, x1, RJ, RJp;
	
	for(int L = 1; L <= N; L++) {

		int n = N-L+1;
		x0 = 0.0;

		for(int it = 1;  it <= IT_MAX; it++) {

				RJ = licz_r(a, b, n, x0);
				RJp = licz_r(b, c, n-1, x0);
				x1 = x0 - (RJ/RJp);
				fprintf(file, "%d | %d | %g | %g | %g |\n", L, it, x1, RJ, RJp);
				if(fabs(x1-x0) < 10e-7) break;
				x0 = x1;
		}

		for(int i = 0; i <= (n-1); i++){
			a[i] = b[i];
		}

		fprintf(file, "\n");

	}

	fclose(file);	
	free(a);
	free(b);
	free(c);
}

//////////////////////////////////////////////////////////////////////////////////