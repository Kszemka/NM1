////////////////////////////////////////////////////////////////////////////////////////
// Method of steepest descent
// using GLS library
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>

#define h 0.1
#define delta 1.0e-4
#define MAX 1000
#define eps1 1.0e-2
#define eps2 1.0e-3


double function(double x, double y);

double df_dx(double x, double y);

double df_dy(double x, double y);


///////////////////////////////////////////////////////////////////////////////////////

int main() {

    FILE* file1 = fopen("eps1.dat", "w");
    FILE* file2 = fopen("eps2.dat", "w");

    double x = -0.75;
    double y = 1.75;

	fprintf(file1, "%f %f\n", x, y);
	fprintf(file2, "%f %f\n", x, y);


    for (int i = 0; i < MAX; i++) {
    	double x1 = x - h * df_dx(x, y);
		double y1 = y - h * df_dy(x, y);
        
		fprintf(file1, "%f %f\n", x1, y1);

	    if (sqrt(pow(x1-x,2)+pow(y1-y,2)) < eps1) break;
	
		x = x1;
		y = y1;
    }

    x = -0.75;
    y = 1.75;

    for (int i = 0; i < MAX; i++) {
        double x1 = x - h * df_dx(x, y);
        double y1 = y - h * df_dy(x, y);

        fprintf(file2, "%f %f\n", x1, y1);

	    if (sqrt(pow(x1-x,2)+pow(y1-y,2)) < eps2) break;

        x = x1;
        y = y1;
    }

    fclose(file1);
    fclose(file2);

}


////////////////////////////////////////////////////////////////////////////////////////////////
//// implementacja funkcji pomocniczych
///////////////////////////////////////////////////////////////////////////////////////////////


double function(double x, double y) {
    return 5.0/2.0 * pow(pow(x,2)-y,2) + pow((1 - x),2);
}

double df_dx(double x, double y) {
    return (function(x + delta, y) - function(x - delta, y)) / (2.0 * delta);
}

double df_dy(double x, double y) {
    return (function(x, y + delta) - function(x, y - delta)) / (2.0 * delta);
}