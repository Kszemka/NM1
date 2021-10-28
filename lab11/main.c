////////////////////////////////////////////////////////////////////////////////////////
// Fast Fourier Transform
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define T 1.0
#define omega 2*M_PI/T

float delta();

void FFT(int x, FILE* file);

//////////////////////////////////////////////////////////////////////////////////////////////

int main() {

    srand(time(NULL));

    FILE* file1 = fopen("k8.dat", "w");
    FILE* file2 = fopen("k10.dat", "w");
    FILE* file3 = fopen("k12.dat", "w");

    FFT(8, file1);
    FFT(10, file2);
    FFT(12, file3);

    fclose(file1);
    fclose(file2);
    fclose(file3);
    
}


////////////////////////////////////////////////////////////////////////////////////////////////

float delta() {
    return (float)rand()/(RAND_MAX+1.0) - 0.5;
}

void FFT(int x, FILE* file){

    const int N = pow(2, x);
    const float dt = (float)(3*T/N);
    const float sigma = T/20.0;

    double* f = malloc(2*N*sizeof(double));
    double* dane = malloc(2*N*sizeof(double));
    double* g1 = malloc(2*N*sizeof(double));
    double* g2 = malloc(2*N*sizeof(double));

   for(int i = 0; i < N; i++) {
        float ti = dt*i;
        f[2*i] = dane[2*i] = sin(omega*ti) +  sin(2*omega*ti) + sin(3*omega*ti) + delta();
        f[2*i + 1] = dane[2*i + 1] = g1[2*i + 1] = g2[2*i + 1] = 0;
        g1[2*i] = g2[2*i] = (1.0 * exp(-1.0*pow(ti,2)/(2*pow(sigma,2)))) / (sigma * sqrt(2*M_PI));
    }
            
    gsl_fft_complex_radix2_forward(f, 1, N);
    gsl_fft_complex_radix2_forward(g1, 1, N);
    gsl_fft_complex_radix2_backward(g2, 1, N);  

    for (int i = 0; i < N; i++) {
        float a1 = f[2*i];
        float b1 = f[2*i+1];
        float a2 = g1[2*i] + g2[2*i];
        float b2 = g1[2*i+1] + g2[2*i+1];
        f[2*i] = a1 * a2 - b1 * b2;
        f[2*i+1] = a1 * b2 + a2 * b1;
    }

    gsl_fft_complex_radix2_backward(f, 1, N);
    double max = fabs(f[1]);

    for(int i = 1; i <= N; i++) {
        if(fabs(f[2*i])>max) max = fabs(f[2*i]);
    }

    for(int i = 0; i <= N; i++) {
        double ti = dt*i;
        fprintf(file, "%f %f\n", ti, dane[2*i]);
    }
     
    fprintf(file, "\n\n");

    for (int i = 0; i <= N; i++) {
        double ti = dt*i;
        fprintf(file, "%f %f\n", ti, f[2*i] * 2.5 / max);
    }
}