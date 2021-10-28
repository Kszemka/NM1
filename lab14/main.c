////////////////////////////////////////////////////////////////////////////////////////
// Generating pseudorandom numbers with specific distribution
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gen_1(long int m);

double gen_2(long int m);

void U_1_2(int N, long int m, double (*generator)(long int), FILE * file1, FILE * file2);

double gen_trian(int delta, int ni);

void triangle(int N, FILE * file);

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(){

	FILE* f1 = fopen("U.dat", "w");
	FILE* f2 = fopen("U_hist.dat", "w");
	FILE* f3 = fopen("T_hist.dat", "w");

	int N = 10000;

	U_1_2(N, pow(2,15), gen_1, f1, f2);
	fprintf(f1, "\n\n");
	fprintf(f2, "\n\n");
	U_1_2(N, pow(2,32), gen_2, f1, f2);

	N = 1000;

	triangle(N, f3);

	fclose(f1);
	fclose(f2);
	fclose(f3);

}

////////////////////////////////////////////////////////////////////////////////////////////////

double gen_1(long int m){
    static long int x = 10;
    int a = 123;
    x = (a*x+1)%m;
    return x/(m+1.0);
}

double gen_2(long int m){
    static long int x = 10;
    int a = 69069;
    x = (a*x+1)%m;
    return x/(m+1.0);
}

void U_1_2(int N, long int m, double (*generator)(long int), FILE *file1, FILE * file2){
    double *x = malloc(N*sizeof(double));
    int K = 12;
    double sum = 0;
    double var = 0;

    for(int i=0; i<N; i++){
        x[i] = generator(m);
        sum+=x[i];
    }

    for(int i = 0; i < N-2; ++i){
        fprintf(file1, "%f\t%f\n", x[i], x[i+1]);
    }

    //srednia
    double ni = (1/N)*sum;
    printf("średnia = %g\n", ni);

    for(int i=0; i<N; i++){
        var+=(x[i]-ni)*(x[i]-ni);
    }

    //odchylenie
    double sd = sqrt((1/N)*var);
    printf("odchylenie = %g\n", sd);


    //hist prawdopodobienstwa
    double x_min = 0.0;
    double x_max = 1.0;
    double dx = (x_max-x_min)/K;

    double *nj = malloc(K*sizeof(double));
    for(int i=0; i<K; i++){
    	nj[(int)((x[i]-x_min)/dx)]++;
    	fprintf(file2, "%f\t%f\n",  (x_min+i*dx+x_min+(i+1)*dx)/2, nj[i]/N);
    	x_min+=dx;
    }

    free(x);
}


double gen_trian(int delta, int ni){
	return ni+(gen_2(pow(2,32))+gen_2(pow(2,32))-1.0)*delta;
}

void triangle(int N, FILE * file){
	int K = 10;

	double *x = malloc(N*sizeof(double));
	for(int i=0; i<N; i++){
		x[i] = gen_trian(3, 4);
	}

	double x_min = 1.0; 
	double x_max = 7.0;
	double sum = 0.0;
    double p_i, x_jmin, x_jmax, f_max, f_min;
	double dx = (x_max-x_min)/K;

    double *nj = malloc(K*sizeof(double));
    for(int i=0; i<K; i++){
    	x_jmin = x_min+i*dx; 
    	x_jmax = x_min+(i+1)*dx;
    	if(x_jmin <= 4) f_min = -1.0/9.0*(-1.0*(x_jmin/2.0)+4*x_jmin)+(x_jmin/3.0);
    	else f_min = -1.0/9.0*((x_jmin/2.0)-4*x_jmin+16.0)+(x_jmin/3.0);
    	
    	if(x_jmax<=4) f_max = -1.0/9.0*(-1.0*(x_jmax/2.0)+4*x_jmax)+(x_jmax/3.0);
    	else f_max = -1.0/9.0*((x_jmax/2.0)-4*x_jmax+16.0)+(x_jmax/3.0);	

    	p_i = f_max - f_min;
    	nj[(int)((x[i]-x_min)/dx)]++;

    	sum+=(nj[i]-N*p_i)*(nj[i]-N*p_i)/(N*p_i);

    	fprintf(file, "%f\t%f\t%f\n",  (x_min+i*dx+x_min+(i+1)*dx)/2, nj[i]/N, p_i);
    	x_min+=dx;
    }

    printf("statystyka testowa: %f\n", sum);
}
