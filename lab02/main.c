////////////////////////////////////////////////////////////////////////////////////////
// Calculating inverse of a matrix and condition number with LU decomposition 
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define N 4
#define delta 2

void fprint_matrix(gsl_matrix * A, FILE * file);

void create_A(gsl_matrix *A);

void create_U(gsl_matrix *U, gsl_matrix *A);

void create_L(gsl_matrix *L, gsl_matrix *A);

float det_matrix(gsl_matrix *A);

void scalar_matrix(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);

float find_max(gsl_matrix *A);


///////////////////////////////////////////////////////////////////////////////////////

int main(){
   
    //matrix A
    gsl_matrix *A = gsl_matrix_calloc(N,N);
    gsl_permutation *p = gsl_permutation_calloc(N);
    
    FILE *file_ptr;
    file_ptr = fopen("out.txt","w");
    fprintf(file_ptr, "Macierz A:\n");

    create_A(A);    
    fprint_matrix(A, file_ptr);

    //matrix LU
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);

    gsl_matrix *U = gsl_matrix_calloc(N,N);
    gsl_matrix *L = gsl_matrix_calloc(N,N);

    create_U(U,A);
    create_L(L,A);

    fprintf(file_ptr, "U macierz: \n");
    fprint_matrix(U, file_ptr);

    fprintf(file_ptr, "L macierz: \n");
    fprint_matrix(L, file_ptr);

    float detA = det_matrix(U)*signum;
    fprintf(file_ptr, "detA: %e\n\n", detA);

    //an invers of A matrix -> A^-1
    gsl_matrix *inverse = gsl_matrix_calloc(N,N);

    gsl_vector *b = gsl_vector_calloc(N);   
    gsl_vector *b1 = gsl_vector_calloc(N);
    gsl_vector *b2 = gsl_vector_calloc(N);
    gsl_vector *b3 = gsl_vector_calloc(N);
    gsl_vector *b4 = gsl_vector_calloc(N);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 0, 1);
    gsl_linalg_LU_solve(A,p,b,b1);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 1, 1);
    gsl_linalg_LU_solve(A,p,b,b2);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 2, 1);
    gsl_linalg_LU_solve(A,p,b,b3);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 3, 1);
    gsl_linalg_LU_solve(A,p,b,b4);


    for (int i = 0; i < N; i++){
        gsl_matrix_set(inverse, 0, i, gsl_vector_get(b1, i));
        gsl_matrix_set(inverse, 1, i, gsl_vector_get(b2, i));
        gsl_matrix_set(inverse, 2, i, gsl_vector_get(b3, i));
        gsl_matrix_set(inverse, 3, i, gsl_vector_get(b4, i));
    }

    fprintf(file_ptr, "Macierz A^(-1): \n");
    fprint_matrix(inverse, file_ptr);

    //matrix AA = A*A^-1
    create_A(A);
    gsl_matrix *AA = gsl_matrix_calloc(N,N);
    scalar_matrix(A,inverse,AA);
    fprintf(file_ptr, "Macierz A*A^(-1): \n");
    fprint_matrix(AA, file_ptr);

    //condition number
    float cond = find_max(A)*find_max(inverse);
    fprintf(file_ptr, "Macierz A max: %f\t Macierz A^(-1) max: %f\t Wskaźnik: %f\n", find_max(A), find_max(inverse), cond);

    gsl_matrix_free(A);
    gsl_matrix_free(U);
    gsl_matrix_free(L);
    gsl_matrix_free(inverse);
    gsl_matrix_free(AA);
    gsl_permutation_free(p);
    gsl_vector_free(b);
    gsl_vector_free(b1);
    gsl_vector_free(b2);
    gsl_vector_free(b3);
    gsl_vector_free(b4);

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

void create_A(gsl_matrix *A){
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            gsl_matrix_set(A,i,j,1.0/(i+j+delta));
        }
    }
}

void create_U(gsl_matrix *U, gsl_matrix *A){
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            if(i <= j){
                gsl_matrix_set(U, i, j, gsl_matrix_get(A,i,j));
            }
            else{
                gsl_matrix_set(U, i, j, 0);
            }
        }
    }
}

void create_L(gsl_matrix *L, gsl_matrix *A){
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            if(j == i){
                gsl_matrix_set(L, i, j, 1);
            }
            else if(i > j){
                gsl_matrix_set(L, i, j, gsl_matrix_get(A,i,j));
            }
            else{
                gsl_matrix_set(L, i, j, 0);

            }
        }   
    }
}

float det_matrix(gsl_matrix *A){
    float det = 1.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            if(i==j){
                det *= gsl_matrix_get(A,i,j);
            }    
        }
    }
    return det;
}

void scalar_matrix(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C){
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            float temp = 0.0;
            for (int k = 0; k < N; k++){
                temp += gsl_matrix_get(A,i,k) * gsl_matrix_get(B,k,j);
            }
            gsl_matrix_set(C,i,j,temp);
        }
    }
}

float find_max(gsl_matrix *A){
    float max = 0.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            if(gsl_matrix_get(A,i,j)>max){
                max = gsl_matrix_get(A,i,j);
            }
        }
    }
    return max;
}
