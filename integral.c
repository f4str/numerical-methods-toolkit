#include "integral.h"

double left_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double delta = (b - a) / n;
	
	for (i = 0; i < n; ++i) {
		sum += f(a + i * delta);
	}
	return delta * sum;
}

double right_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double delta = (b - a) / n;
	
	for (i = 1; i <= n; ++i) {
		sum += f(a + i * delta);
	}
	return delta * sum;
}

double midpoint_rule(double (*f)(double), double a, double b, int n) {
	double i;
	double sum = 0;
	double delta = (b - a) / n;
	
	for (i = 0.5; i < n; ++i) {
		sum += f(a + i * delta);
	}
	return delta * sum;
}

double trapezoidal_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double delta = (b - a) / n;
	
	sum += f(a);
	for (i = 1; i <= n; ++i) {
		sum += 2 * f(a + i * delta);
	}
	sum += f(b);
	
	return delta / 2 * sum;
}

double simpson_rule(double (*f)(double), double a, double b, int n) {
	int i;
	double sum = 0;
	double delta = (b - a) / n;
	
	sum += f(a);
	for (i = 1; i <= n; i += 2) {
		if (i % 2 == 0) {
			sum += 2 * f(a + i * delta);
		}
		else {
			sum += 4 * f(a + i * delta);
		}	
	}
	sum += f(b);
	
	return delta / 3 * sum;
}
