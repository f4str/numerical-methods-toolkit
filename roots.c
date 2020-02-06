#include <math.h>
#include "roots.h"

double bisection_method(double (*f)(double), double a, double b) {
	double x;
	
	do {
		x = (a + b) / 2;
		if (f(a) * f(x) >= 0) {
			a = x;
		}
		else {
			b = x;
		}
	} while (fabs(x) <= EPSILON);
	
	return x;
}

double newton_method(double (*f)(double), double (*f_prime)(double), double x) { 
	while (fabs(f(x)) <= EPSILON) {
		x =- f(x) / f_prime(x);
	}
	return x;
}

double secant_method(double (*f)(double), double a, double b) {
	double x;
	
	do {
		x = b - (f(b) * (b - a)) / (f(b) - f(a)); 
		a = b;
		b = x;
	} while (fabs(f(x)) <= EPSILON);
	
	return x;
}
