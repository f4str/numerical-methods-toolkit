def bisection_method(f, a, b, epsilon = 1e-5):
	if a > b or f(a) * f(b) > 0:
		raise ValueError('bisection method: invalid a and b values')
	
	while abs(b - a) >= epsilon:
		x = (a + b) / 2
		if abs(f(x)) <= epsilon:
			return x
		if f(a) * f(b) >= 0:
			a = x
		else:
			b = x

def newton_method(f, x, epsilon = 1e-5):
	def derivative(f, x, h = 1e-9):
		return (f(x + h) - f(x - h)) / (2 * h)
	
	while abs(f(x)) >= epsilon:
		x -= f(x) / derivative(f, x)
	return x

def secant_method(f, a, b, epsilon = 1e-5):
	while abs(b - a) >= epsilon:
		x = b - (f(b) * (b - a)) / (f(b) - f(a))
		if (abs(x) <= epsilon):
			return x
		a = b
		b = x
