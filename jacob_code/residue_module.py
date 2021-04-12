# -*- coding: utf-8 -*-
import numpy as np
from numpy import sqrt
from itertools import permutations

"""
DISPLAYING RESIDUES
"""

# Displays an input residue matrix inline using Unicode.
def show(arr):
    dim = arr.shape[0],arr.shape[1]
    result = np.zeros(dim,dtype=object)
    for i in range(dim[0]):
        for j in range(dim[1]):
            if arr[i,j] == 0:
                result[i,j] = u'0'
            elif arr[i,j] == 1:
                result[i,j] = u'1'
            elif arr[i,j] == sqrt(2):
                result[i,j] = u'\u221A2'
            else:
                result[i,j] = u'1+\u221A2'
    print(result)
    
"""
CREATING RESIDUE ROWS
"""

# Generates all 14 normalized residue rows up to permutation.
def generate_residue_rows():
    rows = []
    for c1 in [[0,0],[1,0],[0,1],[1,1]]:
        for c2 in [[0,0],[1,0],[0,1],[1,1]]:
            for c3 in [[0,0],[1,0],[0,1],[1,1]]:
                for c4 in [[0,0],[1,0],[0,1],[1,1]]:
                    for c5 in [[0,0],[1,0],[0,1],[1,1]]:
                        for c6 in [[0,0],[1,0],[0,1],[1,1]]:
                            row = [c1[0]+sqrt(2)*c1[1],c2[0]+sqrt(2)*c2[1]
                                   ,c3[0]+sqrt(2)*c3[1],c4[0]+sqrt(2)*c4[1]
                                   ,c5[0]+sqrt(2)*c5[1],c6[0]+sqrt(2)*c6[1]]
                            if normal(row):
                                row.sort(reverse=True)
                                if row not in rows:
                                    rows.append(row)
    rows = np.array(rows[::-1])
    swap(rows,1,6)
    swap(rows,2,4)
    swap(rows,3,4)
    swap(rows,4,10)
    swap(rows,5,8)
    swap(rows,9,10)
    swap(rows,10,12)
    swap(rows,6,8)
    swap(rows,7,9)
    swap(rows,8,11)
    swap(rows,9,12)
    swap(rows,10,12)
    swap(rows,10,11)
    return rows

# Applies normalization constraints to check if a residue array is normal.
def normal(arr):
    sum1,sum2 = 0,0
    for i in range(len(arr)):
        x,y = 0,0
        if (arr[i] == 1) or (arr[i] == 1+sqrt(2)):
            x = 1
        if (arr[i] == sqrt(2)) or (arr[i] == 1+sqrt(2)):
            y = 1
        sum1 += x + 2*y
        sum2 += x*y
    cons1 = sum1 % 4 == 0
    cons2 = sum2 % 2 == 0
    return cons1 & cons2

# Useful function for reordering the rows list
def swap(rows,a,b): 
    temp = rows[a].copy()
    rows[a] = rows[b].copy()
    rows[b] = temp.copy()

"""
PERMUTING RESIDUES
"""

# Takes an array arr as input, returns all permutations of arr.
def gen_perms(arr,dim=None,unique=True):
    if dim == None:
        dim = arr.ndim
    results = []
    if dim == 1:
        perms = np.array(list(permutations(arr)))
        for p in perms:
            if unique:
                if not repeat_array(p,results):
                    results.append(p)
            else:
                results.append(p)
    elif dim == 2:
        row_perms = gen_perms(arr,dim=1,unique=unique)
        for rp in row_perms:
            perm_t = np.transpose(rp)
            col_perms = np.array(list(permutations(perm_t)))
            for cp in col_perms:
                perm = np.transpose(cp)
                if unique:
                    if not repeat_array(perm,results):
                        results.append(perm)
                else:
                    results.append(perm)
    return np.array(results)

# Takes an array new and a matrix previous as input
# Checks if new is contained in previous.
def repeat_array(new,previous):
    for prev in previous:
        if np.allclose(new,prev):
            return True
    return False

"""
RESIDUE UNIQUENESS
"""

# Takes two matrices a and b as input.
# Checks if they are equal up to row and column permutation.
def equal(a,b):
    # Checks if they have the same shape and type
    if (a.shape != b.shape) or (a.dtype != b.dtype):
        return False
    # Checks if they are identical
    if np.allclose(a,b):
        return True
    # Checks if they have the same entries
    a_elems = np.sort(a,axis=None)
    b_elems = np.sort(b,axis=None)
    if not np.allclose(a_elems,b_elems):
        return False
    # Checks if they have the same row sums
    a_row_sums = np.sort(np.sum(a,axis=1))
    b_row_sums = np.sort(np.sum(b,axis=1))
    if not np.allclose(a_row_sums,b_row_sums):
        return False
    # Checks if they have the same column sums
    a_col_sums = np.sort(np.sum(a,axis=0))
    b_col_sums = np.sort(np.sum(b,axis=0))
    if not np.allclose(a_col_sums,b_col_sums):
        return False
    # Checks if they have the same row and column pairs
    if check_pairs(a,b):
        return True
    return False

# Checks if two matrices have the same row and column pairs
def check_pairs(a,b):
    # First check the rows
    row_pairs_match = False
    num_rows = a.shape[0]
    row_pairs_a,row_pairs_b = [],[]
    for i in range(num_rows):
        row_pairs_a.append(pair_str(a[i]))
        row_pairs_b.append(pair_str(b[i]))
    perms_a = [''.join(p) for p in list(permutations(row_pairs_a[0]))]
    for i in range(len(perms_a)):
        if perms_a[i] in row_pairs_b:
            # i corresponds to the potential permutation entry
            temp_b = row_pairs_b.copy()
            is_eq = True
            for j in range(1,len(row_pairs_a)):
                rpa = ''.join(list(permutations(row_pairs_a[j]))[i])
                if rpa not in temp_b:
                    is_eq = False
                    break
            if is_eq:
                row_pairs_match = True
                break
    # Then check the columns
    col_pairs_match = False
    num_cols = a.shape[1]
    col_pairs_a,col_pairs_b = [],[]
    for j in range(num_cols):
        col_pairs_a.append(pair_str(a[:,j]))
        col_pairs_b.append(pair_str(b[:,j]))
    perms_a = [''.join(p) for p in list(permutations(col_pairs_a[0]))]
    for i in range(len(perms_a)):
        if perms_a[i] in col_pairs_b:
            # i corresponds to the potential permutation entry
            temp_b = col_pairs_b.copy()
            temp_b.remove(perms_a[i])
            is_eq = True
            for j in range(1,len(col_pairs_a)):
                rpa = ''.join(list(permutations(col_pairs_a[j]))[i])
                if rpa not in temp_b:
                    is_eq = False
                    break
                else:
                    temp_b.remove(rpa)
            if is_eq:
                col_pairs_match = True
                break
    if row_pairs_match and col_pairs_match:
        return True
    return False

# Generates a string to represent the pairing of an array.
def pair_str(arr,sort=False):
    length = len(arr)
    pair = ''
    for i in range(length):
        if arr[i] == 0:
            pair += '0'
        elif arr[i] == 1:
            pair += '1'
        elif arr[i] == sqrt(2):
            pair += '2'
        elif arr[i] == 1+sqrt(2):
            pair += '3'
    if sort:
        return ''.join(sorted(pair))
    return pair

"""
RESIDUE CONDITIONS
"""

# Takes two residue arrays as input, checks if they are orthogonal.
def orthogonal(top,bot):
    cons1,cons2 = 0,0
    for i in range(len(top)):
        x,y,u,v = 0,0,0,0
        if (top[i] == 1) or (top[i] == 1+sqrt(2)):
            x = 1
        if (top[i] == sqrt(2)) or (top[i] == 1+sqrt(2)):
            y = 1
        if (bot[i] == 1) or (bot[i] == 1+sqrt(2)):
            u = 1
        if (bot[i] == sqrt(2)) or (bot[i] == 1+sqrt(2)):
            v = 1
        cons1 += x*u
        cons2 += x*v + y*u
    cons1 = cons1 % 2 == 0
    cons2 = cons2 % 2 == 0
    if cons1 & cons2:
        return True
    return False


# Takes a matrix mat and array new as input (where the last row of mat is new)
# and checks if the new row is orthogonal to all previous rows in mat.
def new_row_orthogonal(mat,new):
    for i in range(len(mat)-1):
        if not orthogonal(mat[i],new):
            return False
    return True

# Takes a matrix mat as input and checks if every row is orthogonal to every other row.
# Note: The first row is by default excluded since it is assumed that it is orthogonal to all subsequent rows.
def rows_orthogonal(mat,first=False):
    start = 1
    if first:
        start = 0
    for i in range(start,len(mat)-1):
        for j in range(i+1,len(mat)):
            if not orthogonal(mat[i],mat[j]):
                return False
    return True

# Takes a residue array arr as input, returns a final residue entry that normalizes arr.
# Note: If arr cannot be normalized, -1 is returned.
def normalize(arr):
    sum1,sum2 = 0,0
    for i in range(len(arr)):
        x,y = 0,0
        if (arr[i] == 1) or (arr[i] == 1+sqrt(2)):
            x = 1
        if (arr[i] == sqrt(2)) or (arr[i] == 1+sqrt(2)):
            y = 1
        sum1 += x + 2*y
        sum2 += x*y
    for last in [0,1,sqrt(2),1+sqrt(2)]:
        x,y = 0,0
        if (last == 1) or (last == 1+sqrt(2)):
            x = 1
        if (last == sqrt(2)) or (last == 1+sqrt(2)):
            y = 1
        cons1 = (sum1 + (x + 2*y)) % 4 == 0
        cons2 = (sum2 + (x*y)) % 2 == 0
        if cons1 & cons2:
            return last
    return -1

# Takes a residue array arr as input, returns a list of the two final entries that normalize arr.
# Note: If arr cannot be normalized, -1 is returned.
def normalize_two(arr):
    sum1,sum2 = 0,0
    for i in range(len(arr)):
        x,y = 0,0
        if (arr[i] == 1) or (arr[i] == 1+sqrt(2)):
            x = 1
        if (arr[i] == sqrt(2)) or (arr[i] == 1+sqrt(2)):
            y = 1
        sum1 += x + 2*y
        sum2 += x*y
    last = []
    for last1 in [0,1,sqrt(2),1+sqrt(2)]:
        for last2 in [0,1,sqrt(2),1+sqrt(2)]:
            x,y,u,v = 0,0,0,0
            if (last1 == 1) or (last1 == 1+sqrt(2)):
                x = 1
            if (last1 == sqrt(2)) or (last1 == 1+sqrt(2)):
                y = 1
            if (last2 == 1) or (last2 == 1+sqrt(2)):
                u = 1
            if (last2 == sqrt(2)) or (last2 == 1+sqrt(2)):
                v = 1
            cons1 = (sum1 + (x + 2*y) + (u + 2*v)) % 4 == 0
            cons2 = (sum2 + (x*y) + (u*v)) % 2 == 0
            if cons1 & cons2:
                last.append([last1,last2])
    return last

# Takes a residue matrix as input, checks if every column obeys normalization constraints.
def cols_normal(mat):
    for j in range(mat.shape[1]):
        col = mat[:,j]
        if not normal(col):
            return False
    return True

# Takes a residue matrix as input, checks if every column is orthogonal to every other column.
# In this case it cannot be assumed that all columns are orthogonal to the first, so first is set to True.
def cols_orthogonal(mat):
    if rows_orthogonal(np.transpose(mat),first=True):
        return True
    return False

"""
SAVING/READING RESIDUES
"""

# Saves the given array of residues to a txt file.
def save_residues(res_arr, filename):
    res_arr_2d = res_arr.reshape(res_arr.shape[0], -1) 
    np.savetxt(filename,res_arr_2d)

# Reads an array of residues from a given txt file.
def read_residues(filename):
    res_arr_loaded = np.loadtxt(filename)
    res_arr = res_arr_loaded.reshape(res_arr_loaded.shape[0], 6, 6)
    return res_arr
