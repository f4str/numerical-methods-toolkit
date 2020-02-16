#pragma once

static const double EPSILON = 1e-5;

double bisection_method(double (*f)(double), double a, double b);
double newton_method(double (*f)(double), double x);
double secant_method(double (*f)(double), double a, double b);
