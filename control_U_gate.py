# -*- coding: UTF-8 -*-
import cirq
import numpy as np
import math
from gate import Gate

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
#cnot
cnot = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, 1, 0]])          
#Identity matrix
I = np.eye(2)

# ConU_example = np.array([[1, 0, 0, 0],
#                         [0, 1, 0, 0],
#                         [0, 0, 1, 1],
#                         [0, 0, 1, 0]])

'''random generate 4x4 unitary matrix'''
ran_matrix = cirq.testing.random_unitary(4)
print(ran_matrix)

"""check if this gate is Control U gate"""
def check_Con_U(matrix:np.array):

    state00 = np.dot(matrix, v00)
    state01 = np.dot(matrix, v01)
    state10 = np.dot(matrix, v10)
    state11 = np.dot(matrix, v11)

    if((state00==v00).all() and (state01==v01).all() and state10[1]==0 
        and state10[1]==0 and state11[0]==0 and state11[1]==0):
        a = state10[2]
        b = state10[3]
        c = state11[2]
        d = state11[3]
        a_ = a.conjugate()
        b_ = b.conjugate()
        # print(a, b, c, d, a_, b_)
        test1 = a_*c+b_*d
        # print(test1)
        if(test1 == 1):
            print('This gate is a Control U gate.')
    else:
        print('This gate is not a Control U gate.')

    return 

'''Decompose Control U Gate with two CNOT gate and some single qubit gates'''
def decompose_Con_U(a, b, c, d):
    Gc = np.dot(Gate('S',b).to_matrix(),Gate('R',c).to_matrix())
    Ga = Gate('S',(d-b)/2).to_matrix()
    Gb = np.dot(Gate('R', (-c)).to_matrix(),Gate('S',(b-d)/2).to_matrix())
    Ge = Gate('E', a).to_matrix() 
    
    Oper1 = np.kron(I, Gc)
    Oper2 = np.kron(I, Gb)
    Oper3 = np.kron(Ge, Ga)


    return Oper1.dot(cnot).dot(Oper2).dot(cnot).dot(Oper3)

check_Con_U(ran_matrix)
control_U = decompose_Con_U(0.1, 0, 0, 0.1)
print(control_U)


  
