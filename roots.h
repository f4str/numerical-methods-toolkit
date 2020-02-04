#pragma once

double bisection_method(double (*f)(double), double a, double b, double epsilon);
double newton_method(double (*f)(double), double (*f_prime)(double), double x, double epsilon);
double secant_method(double (*f)(double), double a, double b, double epsilon);
