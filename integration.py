def left_rule(f, a, b, n):
	h = (b - a) / n
	area = 0
	for i in range(n):
		area += f(a + i * h)
	return h * area

def right_rule(f, a, b, n):
	h = (b - a) / n
	area = 0
	for i in range(1, n + 1):
		area += f(a + i * h)
	return h * area

def midpoint_rule(f, a, b, n):
	h = (b - a) / n
	area = 0
	for i in range(n):
		area += f(a + (i + 0.5) * h)
	return h * area

def trapezoidal_rule(f, a, b, n):
	h = (b - a) / n
	area = 0
	
	area += f(a)
	for i in range(1, n):
		area += 2 * f(a + i * h)
	area += f(b)
	return h / 2 * area

def simpsons_rule(f, a, b, n):
	h = (b - a) / n
	area = 0
	
	area += f(a)
	for i in range(1, n):
		if i % 2 == 1:
			area += 4 * f(a + i * h)
		else:
			area += 2 * f(a + i * h)
	area += f(b)
	return h / 3 * area
	