////////////////////////////////////////////////////////////////////////////////////////
// Numerical integration
// using internal library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "num_rec/nrutil.h"
#include "num_rec/nrutil.c"
#include "num_rec/gauleg.c"
#include "num_rec/gaulag.c"
#include "num_rec/gammln.c"
#include "num_rec/gauher.c"

float c1(float x);
float c2(float x);
float c3(float x);
float c2_a(float x);

void Gauss_Leg(int p1, int p2, float(*function)(float), float exact, FILE* f, int ver);

void Gauss_Lau(int p1, int p2,float(*function)(float), float exact, FILE* f);

void Gauss_Her(int p1, int p2, float(*function)(float), float exact, FILE* f);

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(){
    FILE* file = fopen("out.dat", "w");

    float c1_dok = M_PI/3.0;
 	float c2_dok = -0.8700577;
	float c3_dok = 2.0/13.0;

    Gauss_Leg(2, 100, c1, c1_dok, file, 0);
    Gauss_Her(2, 100, c2, c2_dok, file);
    Gauss_Leg(2, 100, c2_a, c2_dok, file, 1);
    Gauss_Lau(2, 10, c3, c3_dok, file);

    fclose(file);
}


////////////////////////////////////////////////////////////////////////////////////////////////

float c1(float x){
    return 1.0/(x*sqrt((x*x)-1));
}

float c2(float x){
    return 0.5*log(fabs(x));
}

float c3(float x){
    return sin(2*x)*exp(-2*x);
}

float c2_a(float x){
    return log(x)*exp(-1.0*(x*x));
}    

void Gauss_Leg(int p1, int p2, float(*function)(float), float exact, FILE* f, int ver){
    for (int n = p1; n <= p2; n++){
        float* x = vector(1, n);
        float* w = vector(1, n);

        if (ver) gauleg(0, 5, x, w, n);
        else gauleg(1, 2, x, w, n);
        float sum = 0.0;
        for (int i = 1; i <= n; i++){
            sum += function(x[i]) * w[i];
        }
        fprintf(f, "%d %f\n", n, fabs(exact - sum));
        
        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }

    fprintf(f, "\n\n");	
}

void Gauss_Lau(int p1, int p2,float(*function)(float), float exact, FILE* f){
    for (int n = p1; n <= p2; n++){
        float* x = vector(1, n);
        float* w = vector(1, n);
        
        gaulag(x, w, n, 0);
        float sum = 0.0;
        for (int i = 1; i <= n; i++){
            sum += function(x[i]) * w[i];
        }
        fprintf(f, "%d %f\n", n, fabs(exact - sum));
     
        free_vector(x, 1, n);
    	free_vector(w, 1, n);
    }

    fprintf(f, "\n\n");	
}

void Gauss_Her(int p1, int p2, float(*function)(float), float exact, FILE* f){
    for (int n = p1; n <= p2; n = n+2){
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauher(x, w, n);
        float sum = 0.0;
        for (int i = 1; i <= n; i++){
            sum += function(x[i]) * w[i];
        }

        fprintf(f, "%d %f\n", n, fabs(exact - sum));
        
        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }

    fprintf(f, "\n\n");	
}
