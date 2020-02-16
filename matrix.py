def dot_product(x, y):
	if len(x) != len(y):
		raise ValueError('dot product: invalid vector sizes')
	
	product = 0
	for a, b in zip(x, y):
		product += a * b
	return product

def matrix_addition(A, B):
	if len(A) != len(B) or len(A[0]) != len(B[0]):
		raise ValueError('matrix addition: invalid matrix sizes')
	
	rows = len(A)
	cols = len(A[0])
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			C[i][j] = A[i][j] + B[i][j]
	return C

def matrix_subtraction(A, B):
	if len(A) != len(B) or len(A[0]) != len(B[0]):
		raise ValueError('matrix addition: invalid matrix sizes')
	
	rows = len(A)
	cols = len(A[0])
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			C[i][j] = A[i][j] - B[i][j]
	return C

def naive_matrix_multiplication(A, B):
	if len(A[0]) != len(B):
		raise ValueError('naive matrix multiplication: invalid matrix sizes')
	
	rows = len(A)
	cols = len(B[0])
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			for k in range(len(B)):
				C[i][j] += A[i][k] * B[k][j]
	return C

def matrix_split(A):
	if len(A) % 2 != 0 or len(A[0]) % 2 != 0:
		raise ValueError('split matrix: odd matrix not supported')
	
	rows = len(A)
	cols = len(A[0])
	rows_mid = rows // 2
	cols_mid = cols // 2
	
	A_11 = [[A[i][j] for j in range(cols_mid)] for i in range(rows_mid)]
	A_12 = [[A[i][j] for j in range(cols_mid, cols)] for i in range(rows_mid)]
	A_21 = [[A[i][j] for j in range(cols_mid)] for i in range(rows_mid, rows)]
	A_22 = [[A[i][j] for j in range(cols_mid, cols)] for i in range(rows_mid, rows)]
	
	return A_11, A_12, A_21, A_22
	

def strassen_algorithm(A, B):
	if len(A[0]) != len(B):
		raise ValueError('strassen algorithm: invalid matrix sizes')
	
	rows = len(A)
	cols = len(B[0])
	if rows <= 4 or cols <= 4:
		return naive_matrix_multiplication(A, B)
	
	A_11, A_12, A_21, A_22 = matrix_split(A)
	B_11, B_12, B_21, B_22 = matrix_split(B)
	
	M_1 = strassen_algorithm(matrix_addition(A_11, A_22), matrix_addition(B_11, B_22))
	M_2 = strassen_algorithm(matrix_addition(A_21, A_22), B_11)
	M_3 = strassen_algorithm(A_11, matrix_subtraction(B_12, B_22))
	M_4 = strassen_algorithm(A_22, matrix_subtraction(B_21, B_11))
	M_5 = strassen_algorithm(matrix_addition(A_11, A_12), B_22)
	M_6 = strassen_algorithm(matrix_subtraction(A_21, A_11), matrix_addition(B_11, B_12))
	M_7 = strassen_algorithm(matrix_subtraction(A_11, A_22), matrix_addition(B_21, B_22))
	
	C_11 = matrix_addition(matrix_subtraction(matrix_subtraction(M_1, M_4), M_5), M_7)
	C_12 = matrix_addition(M_3, M_5)
	C_21 = matrix_addition(M_2, M_4)
	C_22 = matrix_addition(matrix_addition(matrix_subtraction(M_1, M_2), M_3), M_6)
	
	C = [[0 for _ in range(rows)] for _ in range(cols)]
	for i in range(rows):
		for j in range(cols):
			C[i][j] = C_11[i][j]
			C[i][j + cols] = C_12[i][j]
			C[i + rows][cols] = C_21[i][j]
			C[i + rows][j + cols] = C_22[i][j]
	
	return C
