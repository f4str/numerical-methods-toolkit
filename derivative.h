#pragma once

static const double DELTA = 1e-9;

double forward_difference(double (*f)(double), double x);
double backward_difference(double (*f)(double), double x);
double central_difference(double (*f)(double), double x);
