import math

def dot_product(x, y):
	if len(x) != len(y):
		raise ValueError('dot product: not equal vector sizes')
	
	product = 0
	for a, b in zip(x, y):
		product += a * b
	return product

def matrix_addition(A, B):
	if len(A) != len(B) or len(A[0]) != len(B[0]):
		raise ValueError('matrix addition: not equal matrix sizes')
	
	rows = len(A)
	cols = len(A[0])
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			C[i][j] = A[i][j] + B[i][j]
	return C

def matrix_subtraction(A, B):
	if len(A) != len(B) or len(A[0]) != len(B[0]):
		raise ValueError('matrix subtraction: not equal matrix sizes')
	
	rows = len(A)
	cols = len(A[0])
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			C[i][j] = A[i][j] - B[i][j]
	return C

def standard_matrix_multiplication(A, B):
	if len(A[0]) != len(B):
		raise ValueError('standard matrix multiplication: invalid matrix sizes')
	
	rows = len(A)
	cols = len(B[0])
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			for k in range(len(B)):
				C[i][j] += A[i][k] * B[k][j]
	return C

def strassen_algorithm_square(A, B):
	if len(A) != len(A[0]) or len(B) != len(B[0]) or len(A) != len(B):
		raise ValueError('strassen algorithm square: not square matrices of equal size')
	
	n = len(A)
	if n <= 4:
		return standard_matrix_multiplication(A, B)
	
	mid = n // 2
	
	A_11 = [[0 for _ in range(mid)] for _ in range(mid)]
	A_12 = [[0 for _ in range(mid)] for _ in range(mid)]
	A_21 = [[0 for _ in range(mid)] for _ in range(mid)]
	A_22 = [[0 for _ in range(mid)] for _ in range(mid)]
	
	B_11 = [[0 for _ in range(mid)] for _ in range(mid)]
	B_12 = [[0 for _ in range(mid)] for _ in range(mid)]
	B_21 = [[0 for _ in range(mid)] for _ in range(mid)]
	B_22 = [[0 for _ in range(mid)] for _ in range(mid)]
	
	for i in range(mid):
		for j in range(mid):
			A_11[i][j] = A[i][j]
			A_12[i][j] = A[i][j + mid]
			A_21[i][j] = A[i + mid][j]
			A_22[i][j] = A[i + mid][j + mid]
			
			B_11[i][j] = B[i][j]
			B_12[i][j] = B[i][j + mid]
			B_21[i][j] = B[i + mid][j]
			B_22[i][j] = B[i + mid][j + mid]
	
	M_1 = strassen_algorithm_square(matrix_addition(A_11, A_22), matrix_addition(B_11, B_22))
	M_2 = strassen_algorithm_square(matrix_addition(A_21, A_22), B_11)
	M_3 = strassen_algorithm_square(A_11, matrix_subtraction(B_12, B_22))
	M_4 = strassen_algorithm_square(A_22, matrix_subtraction(B_21, B_11))
	M_5 = strassen_algorithm_square(matrix_addition(A_11, A_12), B_22)
	M_6 = strassen_algorithm_square(matrix_subtraction(A_21, A_11), matrix_addition(B_11, B_12))
	M_7 = strassen_algorithm_square(matrix_subtraction(A_11, A_22), matrix_addition(B_21, B_22))
	
	C_11 = matrix_addition(matrix_subtraction(matrix_subtraction(M_1, M_4), M_5), M_7)
	C_12 = matrix_addition(M_3, M_5)
	C_21 = matrix_addition(M_2, M_4)
	C_22 = matrix_addition(matrix_addition(matrix_subtraction(M_1, M_2), M_3), M_6)
	
	C = [[0 for _ in range(n)] for _ in range(n)]
	
	for i in range(mid):
		for j in range(mid):
			C[i][j] = C_11[i][j]
			C[i][j + mid] = C_12[i][j]
			C[i + mid][j] = C_21[i][j]
			C[i + mid][j + mid] = C_22[i][j]
	
	return C

def strassen_algorithm(A, B):
	A_rows = len(A)
	A_cols = len(A[0])
	B_rows = len(B)
	B_cols = len(B[0])
	
	if A_cols != B_rows:
		raise ValueError('strassen algorithm: invalid matrix sizes')
	
	if A_rows <= 4 or A_cols <= 4:
		return standard_matrix_multiplication(A, B)
	
	n = 2 ** max(int(math.ceil(math.log2(A_rows))), int(math.ceil(math.log2(A_cols))))
	A_padded = [[0 for _ in range(n)] for _ in range(n)]
	B_padded = [[0 for _ in range(n)] for _ in range(n)]
	
	for i in range(n):
		for j in range(n):
			if i < A_rows and j < A_cols:
				A_padded[i][j] = A[i][j]
			if i < B_rows and j < B_rows:
				B_padded[i][j] = B[i][j]
	
	C_padded = strassen_algorithm_square(A_padded, B_padded)
	C = [[0 for _ in range(A_rows)] for _ in range(B_cols)]
	
	for i in range(A_rows):
		for j in range(B_cols):
			C[i][j] = C_padded[i][j]
	
	return C
	