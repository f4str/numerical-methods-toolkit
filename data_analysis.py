import numpy as np


def polynomial_interpolate(x, y):
	M = [[x_i ** i for i in range(len(x))] for x_i in x]
	a = np.linalg.solve(M, y)
	
	def p(t):
		return sum(a_i * t ** i for i, a_i in enumerate(a))
	
	return p


def linear_curve_fit(x, y):
	x = np.array(x)
	y = np.array(y)
	n = len(x)
	
	x_bar = np.mean(x)
	y_bar = np.mean(y)
	s_xx = np.sum(np.square(x)) - n * x_bar * x_bar
	s_xy = np.sum(x * y) - n * x_bar * y_bar
	beta_1 = s_xy / s_xx
	beta_0 = y_bar - beta_1 * x_bar
	
	def f(t):
		return beta_0 + beta_1 * t
	
	return f


def exp_curve_fit(x, y):
	x = np.array(x)
	y = np.log(y)
	n = len(x)
	
	x_bar = np.mean(x)
	y_bar = np.mean(y)
	s_xx = np.sum(np.square(x)) - n * x_bar * x_bar
	s_xy = np.sum(x * y) - n * x_bar * y_bar
	beta = s_xy / s_xx
	alpha = np.exp(y_bar - beta * x_bar)
	
	def f(t):
		return alpha * np.exp(beta * t)
	
	return f
