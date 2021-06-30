import numpy as np


def gaussian_elimination(A, b):
	if len(A) != len(A[0]) or len(A[0]) != len(b):
		raise ValueError('gaussian elimination: invalid matrix or vector sizes')
	
	if not np.any(np.diag(A)):
		raise ValueError('gaussian elimination: no solution')
	
	A = np.array(A)
	b = np.array(b)
	n = len(A)
	
	for row in range(n - 1):
		for i in range(row + 1, n):
			factor = A[i, row] / A[row, row]
			for j in range(row, n):
				A[i, j] -= factor * A[row, j]
			b[i] -= factor * b[row]
	
	x = np.empty(n)
	x[n - 1] = b[n - 1] / A[n - 1, n - 1]
	for row in reversed(range(n - 2)):
		sums = b[row]
		for j in range(row + 1, n):
			sums -= A[row, j] * x[j]
		x[row] = sums / A[row, row]
	return x


def jacobi_method(A, b, n=10):
	if len(A) != len(A[0]) or len(A) != len(b):
		raise ValueError('jacobi method: invalid matrix and vector sizes')
	
	D = np.diagflat(np.diag(A))
	R = np.array(A) - D
	
	D_inv = np.linalg.inv(D)
	b = np.array(b)
	x = np.zeros(len(A))
	
	for _ in range(n):
		x = np.dot(D_inv, b - np.dot(R, x))
	
	return x


def gauss_seidel_method(A, b, n=10):
	if len(A) != len(A[0]) or len(A) != len(b):
		raise ValueError('gauss seidel method: invalid matrix and vector sizes')
	
	L = np.tril(A, k=0)
	U = np.triu(A, k=1)
	
	L_inv = np.linalg.inv(L)
	b = np.array(b)
	x = np.zeros(len(A))
	
	for _ in range(n):
		x = np.dot(L_inv, b - np.dot(U, x))
	
	return x


def successive_over_relaxation_method(A, b, w, n=10):
	if len(A) != len(A[0]) or len(A) != len(b):
		raise ValueError('successive over relaxation method: invalid matrix and vector sizes')
	
	D = np.diagflat(np.diag(A))
	L = np.tril(A, k=-1)
	U = np.triu(A, k=1)
	
	inv = np.linalg.inv(D + w * L)
	x = np.zeros(len(A))
	wb = w * np.array(b)
	M = w * U + (w + 1) * D
	
	for _ in range(n):
		x = np.dot(inv, wb - np.dot(M, x))
	
	return x
