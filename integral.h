#pragma once

double left_rule(double (*f)(double), double a, double b, int n);
double right_rule(double (*f)(double), double a, double b, int n);
double midpoint_rule(double (*f)(double), double a, double b, int n);
double trapezoidal_rule(double (*f)(double), double a, double b, int n);
double simpson_rule(double (*f)(double), double a, double b, int n);
