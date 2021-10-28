////////////////////////////////////////////////////////////////////////////////////////
// Calculating solution of linear equation with conjugate gradient method
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define N 1000
#define m 5

void fprint_matrix(gsl_matrix * A, FILE * file);

void fprint_vector(gsl_vector * v, FILE * file);

void create_A(gsl_matrix *A);

void create_b(gsl_vector *b);

double get_scalar(gsl_vector *v, gsl_vector *x);

gsl_vector * vector_add(gsl_vector *v, gsl_vector *x);

gsl_vector * vector_sub(gsl_vector *v, gsl_vector *x);

gsl_vector * matrix_vector_multiply(gsl_vector *v, gsl_matrix *A);

gsl_vector * vector_double_multiply(gsl_vector *v, double a);


///////////////////////////////////////////////////////////////////////////////////////

int main(){
   
    FILE *file_ptr;
    file_ptr = fopen("out.dat","w");

    //wyznaczenie macierzy A

    gsl_matrix *A = gsl_matrix_calloc(N,N);
    create_A(A);    

    //wyznaczenie wektora wyrazow wolnych

    gsl_vector *b = gsl_vector_calloc(N);
    create_b(b);    

    //zaincjalizowanie wektora startowego

    gsl_vector *x = gsl_vector_calloc(N);
    gsl_vector_set_zero(x);

    //wektor reszt

    gsl_vector *r = gsl_vector_calloc(N);
    gsl_vector_memcpy(r,b);
    gsl_vector_sub(r, matrix_vector_multiply(x, A));

    gsl_vector *v = gsl_vector_calloc(N);
    gsl_vector_memcpy(v,r);


	int k = 0;
	double alfa = 0, beta = 0; 

   	while (sqrt(get_scalar(r,r)) > pow (10, -6)){
   		double tmp_r_scalar = get_scalar(r,r);
   		alfa = tmp_r_scalar/get_scalar(v,matrix_vector_multiply(v,A));
   		vector_add(x, vector_double_multiply(v,alfa));
   		vector_sub(r, vector_double_multiply(matrix_vector_multiply(v, A),alfa));
   		beta = get_scalar(r,r)/tmp_r_scalar;
   		gsl_vector_memcpy(v,vector_double_multiply(v,beta));
   		vector_add(v, r);

    	fprintf(file_ptr, "%d \t %g \t %g \t %g \t %g \n", k, sqrt(get_scalar(r,r)), alfa, beta, sqrt(get_scalar(x,x)));
    	k++;
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(v);
    gsl_vector_free(x);
    gsl_vector_free(r);

    fclose(file_ptr);
 }


////////////////////////////////////////////////////////////////////////////////////////////////


 void fprint_matrix(gsl_matrix * A, FILE * file){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            fprintf(file, "%f ", gsl_matrix_get(A,i,j));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

 void fprint_vector(gsl_vector * v, FILE * file){
    for(int i = 0; i < N; i++){
        fprintf(file, "%f \n", gsl_vector_get(v,i));
    }

}

void create_A(gsl_matrix *A){
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            if (abs(i-j)>m) gsl_matrix_set(A,i,j,0);
            else gsl_matrix_set(A,i,j,1.0/(1+abs(i-j)));
        }
    }
}

void create_b(gsl_vector *b){
    for(int i = 0; i < N; i++) {
       gsl_vector_set(b,i,i+1);
    }
}

double get_scalar(gsl_vector *v, gsl_vector *x){
    double temp = 0;
    for(int i = 0; i < N; i++){
        temp += gsl_vector_get(v,i) * gsl_vector_get(x,i);
    }
    return temp;
}

gsl_vector * vector_add(gsl_vector *v, gsl_vector *x){	
	for (int i = 0; i < N; i++){
		gsl_vector_set(v, i, gsl_vector_get(v,i)+gsl_vector_get(x,i));		
	}
	return v;
}

gsl_vector * vector_sub(gsl_vector *v, gsl_vector *x){
	for (int i = 0; i < N; i++){
		gsl_vector_set(v, i, gsl_vector_get(v,i)-gsl_vector_get(x,i));		
	}
	return v;
}

gsl_vector * matrix_vector_multiply(gsl_vector *v, gsl_matrix *A){
	gsl_vector *out = gsl_vector_calloc(N);
	gsl_vector_set_zero(out);

	for(int i = 0; i < N; i++){
		double temp = 0;
		for(int j = 0; j<N; j++){
			temp += gsl_matrix_get(A,i,j)* gsl_vector_get(v,j);
		}
		gsl_vector_set(out, i, temp);
	}
	return out;
}

gsl_vector * vector_double_multiply(gsl_vector *v, double a){
	for (int i = 0; i < N; i++){
		gsl_vector_set(v, i, gsl_vector_get(v,i)*a);
	}
	return v;
}
