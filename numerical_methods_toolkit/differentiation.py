def forward_difference(f, x, h=1e-9):
    return (f(x + h) - f(x)) / h


def backward_difference(f, x, h=1e-9):
    return (f(x) - f(x - h)) / h


def central_difference(f, x, h=1e-9):
    return (f(x + h) - f(x - h)) / (2 * h)
