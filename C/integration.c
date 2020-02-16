#include "integration.h"

double left_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double h = (b - a) / n;
	
	for (i = 0; i < n; ++i) {
		sum += f(a + i * h);
	}
	return h * sum;
}

double right_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double h = (b - a) / n;
	
	for (i = 1; i <= n; ++i) {
		sum += f(a + i * h);
	}
	return h * sum;
}

double midpoint_rule(double (*f)(double), double a, double b, int n) {
	double i;
	double sum = 0;
	double h = (b - a) / n;
	
	for (i = 0.5; i < n; ++i) {
		sum += f(a + i * h);
	}
	return h * sum;
}

double trapezoidal_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double h = (b - a) / n;
	
	sum += f(a);
	for (i = 1; i < n; ++i) {
		sum += 2 * f(a + i * h);
	}
	sum += f(b);
	
	return h / 2 * sum;
}

double simpson_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double h = (b - a) / n;
	
	sum += f(a);
	for (i = 1; i < n; ++i) {
		if (i % 2 == 0) {
			sum += 2 * f(a + i * h);
		}
		else {
			sum += 4 * f(a + i * h);
		}	
	}
	sum += f(b);
	
	return h / 3 * sum;
}
