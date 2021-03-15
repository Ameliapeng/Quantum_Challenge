import cirq
import numpy as np

#statesvector
v00 = np.array([[1],
                [0],
                [0],
                [0]])
v01 = np.array([[0],
                [1],
                [0],
                [0]])
v10 = np.array([[0],
                [0],
                [1],
                [0]])  
v11 = np.array([[0],
                [0],
                [0],
                [1]])              
#CNOT
cnot1 = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, 1, 0]])
cnot2 = np.array([[1, 0, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, 1, 0],
                  [0, 1, 0, 0]])
#SWAP
swap_standard = np.array([[1, 0, 0, 0],
                        [0, 0, 1, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1]])

'''random generate 4x4 unitary matrix'''

ran_matrix = cirq.testing.random_unitary(4)
# print(random_uni_matrix)

'''check if this gate is swap gate'''

state00 = np.dot(ran_matrix, v00)
state01 = np.dot(ran_matrix, v01)
state10 = np.dot(ran_matrix, v10)
state11 = np.dot(ran_matrix, v11)

if ((state00 ==v00).all() and (state01 == v10).all() and (state10 == v01).all() and (state11 ==v11).all()):
    print('This gate is a swap gate:\n')
    swap_gate = ran_matrix
    print(swap_gate)


else:
    print('This gate is not a swap gate.\n')

# state00 = np.dot(swap_standard, v00)
# state01 = np.dot(swap_standard, v01)
# state10 = np.dot(swap_standard, v10)
# state11 = np.dot(swap_standard, v11)

# if ((state00 ==v00).all() and (state01 == v10).all() and (state10 == v01).all() and (state11 ==v11).all()):
#     print('This gate is a swap gate:\n')
#     swap_gate = ran_matrix
#     print(swap_gate)

# else:
#     print('This gate is not a swap gate.\n')

'''Decompose swap gate!'''
swap_decomposition1 = cnot1.dot(cnot2).dot(cnot1)
swap_decomposition2 = cnot2.dot(cnot1).dot(cnot2)




