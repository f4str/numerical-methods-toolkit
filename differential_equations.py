def forward_euler(f, t0, y0, n, h=0.01):
	t = [t0]
	y = [y0]
	t_i = t0
	y_i = y0
	
	for _ in range(n):
		t_i += h
		y_i += h * f(t_i, y_i)
		
		t.append(t_i)
		y.append(y_i)
	
	return t, y


def heun_method(f, t0, y0, n, h=0.01):
	t = [t0]
	y = [y0]
	t_i = t0
	y_i = y0
	
	for _ in range(n):
		t_i += h
		y_pred = y_i + h * f(t_i, y_i)
		y_i += h / 2 * (f(t_i, y_i) + f(t_i, y_pred))
		
		t.append(t_i)
		y.append(y_i)
	
	return t, y


def runge_kutta_method(f, t0, y0, n, h=0.01):
	t = [t0]
	y = [y0]
	t_i = t0
	y_i = y0
	
	for _ in range(n):
		t_i += h
		k1 = f(t_i, y_i)
		k2 = f(t_i + h / 2, y_i + h / 2 * k1)
		k3 = f(t_i + h / 2, y_i + h / 2 * k2)
		k4 = f(t_i + h, y_i + h * k3)
		y_i += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
		
		t.append(t_i)
		y.append(y_i)
	
	return t, y
