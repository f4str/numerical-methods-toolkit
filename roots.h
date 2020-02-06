#pragma once

static const double EPSILON = 0.0001;

double bisection_method(double (*f)(double), double a, double b);
double newton_method(double (*f)(double), double (*f_prime)(double), double x);
double secant_method(double (*f)(double), double a, double b);
