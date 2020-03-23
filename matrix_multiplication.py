import math


def dot_product(x, y):
	if len(x) != len(y):
		raise ValueError('dot product: not equal vector sizes')
	return sum(x_i * y_i for x_i, y_i in zip(x, y))


def matrix_addition(A, B):
	if len(A) != len(B) or len(A[0]) != len(B[0]):
		raise ValueError('matrix addition: not equal matrix sizes')
	
	return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def matrix_subtraction(A, B):
	if len(A) != len(B) or len(A[0]) != len(B[0]):
		raise ValueError('matrix subtraction: not equal matrix sizes')
	
	return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def matrix_vector_multiplication(A, x):
	if len(A[0]) != len(x):
		raise ValueError('matrix vector product: invalid matrix or vector sizes')
	
	return [sum(a * x_i for a, x_i in zip(row, x)) for row in A]


def standard_matrix_multiplication(A, B):
	if len(A[0]) != len(B):
		raise ValueError('standard matrix multiplication: invalid matrix sizes')
	
	return [[sum(a * b for a, b in zip(row, col)) for col in zip(*B)] for row in A]


def strassen_algorithm_square(A, B, levels):
	if len(A) != len(A[0]) or len(B) != len(B[0]) or len(A) != len(B):
		raise ValueError('strassen algorithm square: not square matrices of equal size')
	
	n = len(A)
	if n <= 4 or levels <= 0:
		return standard_matrix_multiplication(A, B)
	
	mid = n // 2
	
	A_11 = [A[i][:mid] for i in range(mid)]
	A_12 = [A[i][mid:] for i in range(mid)]
	A_21 = [A[i][:mid] for i in range(mid, n)]
	A_22 = [A[i][mid:] for i in range(mid, n)]
	
	B_11 = [B[i][:mid] for i in range(mid)]
	B_12 = [B[i][mid:] for i in range(mid)]
	B_21 = [B[i][:mid] for i in range(mid, n)]
	B_22 = [B[i][mid:] for i in range(mid, n)]
	
	M_1 = strassen_algorithm_square(matrix_addition(A_11, A_22), matrix_addition(B_11, B_22), levels - 1)
	M_2 = strassen_algorithm_square(matrix_addition(A_21, A_22), B_11, levels - 1)
	M_3 = strassen_algorithm_square(A_11, matrix_subtraction(B_12, B_22), levels - 1)
	M_4 = strassen_algorithm_square(A_22, matrix_subtraction(B_21, B_11), levels - 1)
	M_5 = strassen_algorithm_square(matrix_addition(A_11, A_12), B_22, levels - 1)
	M_6 = strassen_algorithm_square(matrix_subtraction(A_21, A_11), matrix_addition(B_11, B_12), levels - 1)
	M_7 = strassen_algorithm_square(matrix_subtraction(A_12, A_22), matrix_addition(B_21, B_22), levels - 1)
	
	C_11 = matrix_addition(matrix_subtraction(matrix_addition(M_1, M_4), M_5), M_7)
	C_12 = matrix_addition(M_3, M_5)
	C_21 = matrix_addition(M_2, M_4)
	C_22 = matrix_addition(matrix_addition(matrix_subtraction(M_1, M_2), M_3), M_6)
	
	C = [c11 + c12 for c11, c12 in zip(C_11, C_12)] + [c21 + c22 for c21, c22 in zip(C_21, C_22)]
	
	return C


def strassen_algorithm(A, B, levels):
	A_rows = len(A)
	A_cols = len(A[0])
	B_rows = len(B)
	B_cols = len(B[0])
	
	if A_cols != B_rows:
		raise ValueError('strassen algorithm: invalid matrix sizes')
	
	if A_rows <= 4 or A_cols <= 4:
		return standard_matrix_multiplication(A, B)
	
	n = 2 ** max(int(math.ceil(math.log2(A_rows))), int(math.ceil(math.log2(A_cols))))
	A_padded = [[A[i][j] if i < A_rows and j < A_cols else 0 for j in range(n)] for i in range(n)]
	B_padded = [[B[i][j] if i < B_rows and j < B_cols else 0 for j in range(n)] for i in range(n)]
	
	C_padded = strassen_algorithm_square(A_padded, B_padded, levels)
	C = [[C_padded[i][j] for j in range(A_rows)] for i in range(B_cols)]
	
	return C
