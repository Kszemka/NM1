////////////////////////////////////////////////////////////////////////////////////////
// Defining 1D wave equation with eigenfunction 
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define L 10
#define n 200
#define N 1

void create_A(gsl_matrix * A, const double delta);

void create_B(gsl_matrix * B, int alfa, const double delta);

void vector_fprint(FILE * file, gsl_vector * eval, int alfa);

void matrix_fprint(FILE * file, gsl_matrix * evec, int alfa, const double delta);

double x_i (int i, const double delta){ return -0.5 * L + (delta * (i+1));}

double ro (double x, int alfa){ return 1 + (4 * alfa * pow(x,2));}


///////////////////////////////////////////////////////////////////////////////////////

int main(){
   
    FILE *file_ptr1;
    file_ptr1 = fopen("eval.dat","w"); 

    FILE *file_ptr2;
    file_ptr2 = fopen("evec.dat","w");

    const double delta = L/(n+1);

    //wyznaczenie macierzy A i B

    gsl_matrix *A = gsl_matrix_alloc(n,n);
    gsl_matrix *B = gsl_matrix_alloc(n,n); 


    //zespolony wektor wartości własnych eval, zespolona macierz evec, wektor pomocniczy w

    gsl_vector * eval = gsl_vector_alloc(n);
    gsl_matrix * evec = gsl_matrix_alloc(n, n);
    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(n);


	for (int alfa = 0; alfa <= 100; alfa += 2){

    	create_A(A, delta); 
		create_B(B, alfa, delta);
		gsl_eigen_gensymmv(A, B, eval, evec, w);
		gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		vector_fprint(file_ptr1, eval, alfa);
		
		if( alfa == 0 || alfa == 100) matrix_fprint(file_ptr2, evec, alfa, delta);

	}

	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_eigen_gensymmv_free(w);

    fclose(file_ptr1);
    fclose(file_ptr2);
 }


////////////////////////////////////////////////////////////////////////////////////////////////


void create_A(gsl_matrix *A, const double delta){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if (i==j) gsl_matrix_set(A,i,j, 2. / pow(delta,2));
            if (i == j+1) gsl_matrix_set(A,i,j,-1. / pow(delta,2));
            if (i == j-1) gsl_matrix_set(A,i,j,-1. / pow(delta,2));
            else gsl_matrix_set(A,i,j, 0);
		}
    }
}

void create_B(gsl_matrix *B, int alfa, const double delta){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
        	double tmp = 0;
            if (i==j) gsl_matrix_set(B,i,j,ro(x_i(i, delta), alfa)/N);
            else gsl_matrix_set(B,i,j,0);
		}
    }
}

void vector_fprint(FILE * file, gsl_vector * eval, int alfa){
	for (int i = 0; i < 6; i++){
		fprintf(file, "%d\t%f\n", alfa, gsl_vector_get(eval, i));			
	}
}

void matrix_fprint(FILE * file, gsl_matrix * evec, int alfa, const double delta){
	for (int i = 0; i < 6; i++){
		fprintf(file, "%f\t", x_i(i, delta));
		for (int j = 0; j < n; j++){
			fprintf(file, "%f\t", gsl_matrix_get(evec, i, j));
		}
		fprintf(file, "\n");		
	}
	fprintf(file, "\n\n");
}