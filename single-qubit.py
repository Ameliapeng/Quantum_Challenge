import cirq

'''random generate 4x4 unitary matrix'''

random_uni_matrix = cirq.testing.random_unitary(4)
# print(random_uni_matrix)

'''decompose 4x4 unitary matrix to two 2x2 unitary matrix(kronecker decompose)'''
decompose = cirq.linalg.kron_factor_4x4_to_2x2s(random_uni_matrix)
factor = decompose[0]
matrix1 = decompose[1]
matrix2 = decompose[2]
print('single-qubit qubits:\n')
print(matrix1)
print(matrix2)

