"""
This module defines the Operator class as well as various useful functions
involving the Operator class.
"""

import numpy as np
from copy import deepcopy

"""
CREATING OPERATORS
"""

# Defines a class for the SO(6) operator
class Operator:
    def __init__(self, k, mat, tc, res, name):
        # Int representing the denominator exponent
        self.k = k
        # (2,6,6) int numpy array for the matrix
        self.matrix = mat
        # Int representing the T count
        self.t_count = tc
        # (2,6,6) int numpy array for the residue
        self.residue = res
        # Str that lists all combined T operators
        self.name = name
        # Str that is column permutation invariant, used for comparison
        self.col_str = col_str(mat)
        
    def __eq__(self,other):
        if (isinstance(other,Operator)):
            return (self.k == other.k and np.all(self.matrix == other.matrix)
                    and self.t_count == other.t_count
                    and np.all(self.residue == other.residue)
                    and self.name == other.name)
        return False

# The col_str is a string used to compare operators
# Checks if they are identical up to column permutation
def col_str(mat):
    arr = []
    for j in range(6):
        string = ''
        for i in range(6):
            string += str(mat[0][i,j]) + ':' + str(mat[1][i,j])
            string += ','
        string = string[:-1]
        index = 0
        for k in range(j):
            if string > arr[k]:
                index += 1
        arr = arr[:index] + [string] + arr[index:]
    result = str(arr)
    return result

"""
DISPLAYING OPERATORS
"""

# Displays the matrix for an SO(6) operator
def show(op):
    result = show_unicode(op.matrix)
    print(f'k = {op.k}')
    print(result)

# Displays the residue for an SO(6) operator
def show_residue(op):
    result = show_unicode(op.residue)
    print(result)
        
# Displays an input matrix inline using Unicode.
def show_unicode(arr):
    result = np.zeros((6,6),dtype=object)
    for i in range(6):
        for j in range(6):
            result[i,j] = ''
            if not ((arr[0][i,j] == 0) and (arr[1][i,j] != 0)):
                result[i,j] += f'{arr[0][i,j]}'
                if arr[1][i,j] > 0:
                    result[i,j] += '+'
            if arr[1][i,j] != 0:
                if abs(arr[1][i,j]) != 1:
                    result[i,j] += f'{arr[1][i,j]}' + '\u221A2'
                elif arr[1][i,j] == 1:
                    result[i,j] += '\u221A2'
                elif arr[1][i,j] == -1:
                    result[i,j] += '-\u221A2'
    return result

"""
MULTIPLYING OPERATORS
"""

# Returns the matrix product of two input operators.
def dot(op1,op2):
    name = op1.name + ' ' + op2.name
    tc = op1.t_count + op2.t_count
    k = op1.k + op2.k
    mat0 = np.matmul(op1.matrix[0],op2.matrix[0]) + 2*np.matmul(op1.matrix[1],op2.matrix[1])
    mat1 = np.matmul(op1.matrix[0],op2.matrix[1]) + np.matmul(op1.matrix[1],op2.matrix[0])
    mat = np.array([mat0,mat1],dtype=int)
    pat = get_res(mat)
    result = Operator(k,mat,tc,pat,name)
    result = k_reduce(result)
    return result

# Returns the residue of an SO(6) operator that is taken as input
def get_res(mat):
    return np.array([mat[0] % 2,mat[1] % 2],dtype=int)

# Returns a copy of the input operator that has been k-reduced until it is 
# no longer divisible by sqrt(2).
def k_reduce(op):
    op_copy = deepcopy(op)
    while np.all(op_copy.matrix[0] % 2 == 0):
        op_copy = k_reduce_once(op_copy)
    return op_copy

# Returns a copy of the input operator that has been divided by sqrt(2), 
# and k reduced by 1.
def k_reduce_once(op):
    name = op.name
    tc = op.t_count
    k = op.k - 1
    mat0 = np.copy(op.matrix[1])
    mat1 = np.copy(op.matrix[0])/2
    mat = np.array([mat0,mat1],dtype=int)
    pat = get_res(mat)
    return Operator(k,mat,tc,pat,name)

# NOTE: This turned out to be unnecessary, so dot2() was renamed to just dot()
# # Returns the matrix product of an arbitrary number of input operators, 
# # multiplied from left to right.
# def dot(*operators):
#     result = deepcopy(operators[0])
#     for op in operators[1:]:
#         result = dot2(result,op)
#     return result

"""
COMPARING OPERATORS
"""

# Checks if two operators have SO(6) operators that are identical up to column 
# (not row) permutation.
def is_perm(op1,op2):
    # Checks if they have the same denominator exponent
    if op1.k != op2.k:
        return False
    # Checks if they have the same column strings (which are column-permutation invariant)
    if op1.col_str == op2.col_str:
        return True
    return False

# NOTE: This method is comparatively slow, so it does not need to be used.
# Checks if two operators have SO(6) operators that are identical 
# up to column (not row) permutation via matrix multiplication.
# def is_perm_alt(op1,op2):
#     op1_mat = op1.matrix
#     op1_mat_trans = np.array([np.transpose(op1_mat[0]),np.transpose(op1_mat[1])])
#     op1_trans = operator(op1.k,op1_mat_trans,op1.t_count,get_res(op1_mat_trans),op1.name)
#     prod = dot(op1_trans,op2)
#     if prod.k == 0:
#         return True
#     return False

"""
SAVING/READING OPERATORS
"""

# Saves the given array of operators to a txt file.
def save_operators(op_arr, filename):
    with open(filename, 'x') as writer:
        writer.write(str(len(op_arr)) + '\n')
        for op in op_arr:
            writer.write(str(op.k) + '\n')
            mat_str = str(np.reshape(op.matrix,72).tolist())[1:-1]
            writer.write(mat_str + '\n')
            writer.write(str(op.t_count) + '\n')
            pat_str = str(np.reshape(op.pattern,72).tolist())[1:-1]
            writer.write(pat_str + '\n')
            writer.write(op.name + '\n')

# Reads an array of operators from a given txt file.
def read_operators(filename):
    with open(filename, 'r') as reader:
        result = []
        length = int(reader.readline()[:-1])
        for index in range(length):
            k = int(reader.readline()[:-1])
            mat_str = reader.readline()[:-1]
            mat = np.reshape(np.fromstring(mat_str,dtype=int,sep=', '), (2,6,6))
            tc = int(reader.readline()[:-1])
            pat_str = reader.readline()[:-1]
            pat = np.reshape(np.fromstring(pat_str,dtype=int,sep=', '), (2,6,6))
            name = reader.readline()[:-1]
            op = Operator(k, mat, tc, pat, name)
            result.append(op)
        result = np.array(result)
        return result
