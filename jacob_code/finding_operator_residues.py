"""
This program defines a method for determining the number of unique residues
contained in a given array of operators.
"""

import numpy as np
from numpy import sqrt
from time import time
from copy import deepcopy
import operator_module as opm
import residue_module as rem

# The sum_sort() function is an attempt at sorting residues such that they are
# invariant under column/row permutation. Need to find a way to account for 
# rows/columns having the same sum (which is common for residues).
def sum_sort(matrix):
    mat = deepcopy(matrix)
    num_rows = mat.shape[0]
    row_sums = np.zeros(num_rows)
    for i in range(num_rows):
        row_sums[i] = np.sum(mat[i])
    mat = mat[row_sums.argsort()]
    num_cols = mat.shape[1]
    col_sums = np.zeros(num_cols)
    for i in range(num_cols):
        col_sums[i] = np.sum(mat[:,i])
    mat = mat[:,col_sums.argsort()]
    return mat

# Takes an operator residue as input (in the form (2,6,6) dtype=int) 
# and converts it to a residue (in the form (6,6) dtype=float)
def operator_res(pat):
    res = np.array(pat[0]+pat[1]*sqrt(2),dtype=float)
    return res

# Takes an array of operators and returns the number of unique residues that
# are contained in the complete list of residues.
def num_res(Tn,name):
    start = time()
    # res_dict has a key associated with the index of each entry in residues,
    # and the value is the number of operators that have that residue
    res_dict = dict.fromkeys(range(len(residues)),0)
    for op in Tn:
        if op.k == 1:
            continue
        op_res = operator_res(op.residue)
        for i in range(len(residues)):
            res = residues[i]
            if rem.equal(res,op_res):
                res_dict[i] += 1
                break
    num_res = 0
    for i in range(len(residues)):
        count = res_dict[i]
        if count != 0:
            num_res += 1
    print(f'Residues in {name}: {num_res}')
    print(f'Time taken: {time()-start}')
    return num_res

# Imports the residues array as well as all the arrays of T operators.
residues = rem.read_residues('residues.txt')
T1 = opm.read_operators('T1.txt')
T2 = opm.read_operators('T2.txt')
T3 = opm.read_operators('T3.txt')
# T4 = opm.read_operators('T4.txt')
# T5 = opm.read_operators('T5.txt')

T1_res = num_res(T1,'T1')
T2_res = num_res(T2,'T2')
T3_res = num_res(T3,'T3')

