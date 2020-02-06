#include "derivative.h"

double forward_difference(double (*f)(double), double x) {
	return (f(x + DELTA) - f(x)) / DELTA;
}

double backward_difference(double (*f)(double), double x) {
	return (f(x) - f(x - DELTA)) / DELTA;
}

double central_difference(double (*f)(double), double x) {
	return (f(x + DELTA) - f(x -DELTA)) / (2 * DELTA);
}
