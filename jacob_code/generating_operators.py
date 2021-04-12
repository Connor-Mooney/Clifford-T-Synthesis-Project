# -*- coding: utf-8 -*-
import numpy as np
from time import time
import operator_module as opm

# Generates a list of the T1 operators
def generate_T1(save=False, load=False):
    if load:
        return opm.read_operators('T1.txt')
    T1 = []
    for i in range(6):
        for j in range(i+1,6):
            if (i,j) in ((0,1),(1,2),(3,4),(4,5),(1,0),(2,1),(4,3),(5,4)):
                sign = 1
            else:
                sign = -1
            name = f'T[{i+1},{j+1}]'
            tc = 1
            k = 1
            mat = np.array([np.zeros((6,6)),np.identity(6)],dtype=int)
            mat[0][i,i] = 1
            mat[0][i,j] = -sign
            mat[0][j,i] = sign
            mat[0][j,j] = 1
            mat[1][i,i] = 0
            mat[1][j,j] = 0
            pat = opm.get_res(mat)
            op = opm.Operator(k,mat,tc,pat,name)
            T1.append(op)
    T1 = np.array(T1)
    if save:
        opm.save_operators(T1, 'T1.txt')
    print(f'T1: {len(T1)}')
    return T1

# Generates a list of T2 operators
def generate_T2(save=False, load=False):
    if load:
        return opm.read_operators('T2.txt')
    start = time()
    T2 = []
    for t in T1:
        for op in T1:
            prod = opm.dot(t,op)
            if prod.k == 0:
                continue
            T2.append(prod)
            for prev in np.concatenate((T1,T2[:-1])):
                if opm.is_perm(prev,prod):
                    T2.pop()
    T2 = np.array(T2)
    if save:
        opm.save_operators(T2, 'T2.txt')
    print(f'T2: {len(T2)}')
    print(f'Time taken: {time()-start}')
    k0 = []
    k1 = []
    k2 = []
    for op in T2:
        k = op.k
        if k == 0:
            k0.append(op)
        elif k == 1:
            k1.append(op)
        elif k == 2:
            k2.append(op)
    print(f'k0 = {len(k0)}')
    print(f'k1 = {len(k1)}')
    print(f'k2 = {len(k2)}')
    print(f'total = {len(k0)+len(k1)+len(k2)}')
    return T2

# Generates a list of T3 operators
def generate_T3(save=False, load=False):
    if load:
        return opm.read_operators('T3.txt')
    start = time()
    T3 = []
    for t in T1:
        for op in T2:
            prod = opm.dot(t,op)
            if prod.k == 0:
                continue
            T3.append(prod)
            for prev in T1:
                if opm.is_perm(prev,prod):
                    T3.pop()
                    continue
            for prev in T2:
                if opm.is_perm(prev,prod):
                    T3.pop()
                    continue
            for prev in T3[:-1]:
                if opm.is_perm(prev,prod):
                    T3.pop()
                    continue
    T3 = np.array(T3)
    if save:
        opm.save_operators(T3, 'T3.txt')
    print(f'T3: {len(T3)}')
    print(f'Time taken: {time()-start}')
    k0 = []
    k1 = []
    k2 = []
    k3 = []
    for op in T3:
        k = op.k
        if k == 0:
            k0.append(op)
        elif k == 1:
            k1.append(op)
        elif k == 2:
            k2.append(op)
        elif k == 3:
            k3.append(op)
    print(f'k0 = {len(k0)}')
    print(f'k1 = {len(k1)}')
    print(f'k2 = {len(k2)}')
    print(f'k3 = {len(k3)}')
    print(f'total = {len(k0)+len(k1)+len(k2)+len(k3)}')
    return T3

# Generates a list of T4 operators
def generate_T4(save=False, load=False):
    if load:
        return opm.read_operators('T4.txt')
    start = time()
    T4 = []
    for t in T1:
        for op in T3:
            prod = opm.dot(t,op)
            if prod.k == 0:
                continue
            T4.append(prod)
            for prev in np.concatenate((T1,T2,T3,T4[:-1])):
                if opm.is_perm(prev,prod):
                    T4.pop()
    T4 = np.array(T4)
    if save:
        opm.save_operators(T4, 'T4.txt')
    print(f'T4: {len(T4)}')
    print(f'Time taken: {time()-start}')
    k0 = []
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    for op in T4:
        k = op.k
        if k == 0:
            k0.append(op)
        elif k == 1:
            k1.append(op)
        elif k == 2:
            k2.append(op)
        elif k == 3:
            k3.append(op)
        elif k == 4:
            k4.append(op)
    print(f'k0 = {len(k0)}')
    print(f'k1 = {len(k1)}')
    print(f'k2 = {len(k2)}')
    print(f'k3 = {len(k3)}')
    print(f'k4 = {len(k4)}')
    print(f'total = {len(k0)+len(k1)+len(k2)+len(k3)+len(k4)}')
    return T4

def generate_T5(save=False, load=False):
    if load:
        return opm.read_operators('T5.txt')
    start = time()
    T5 = []
    index = 0
    total = len(T1)*len(T4)
    for t in T1:
        for op in T4:
            index += 1
            if index % int(total * 0.1) == 0:
                dt = time()-start
                print(f'{dt:.2f} {index}')
            prod = opm.dot(t,op)
            if prod.k == 0:
                continue
            T5.append(prod)
            for prev in T5[:-1]:
                if opm.is_perm(prev,prod):
                    T5.pop()
    T5 = np.array(T5)
    if save:
        opm.save_operators(T5, 'T5.txt')
    print(f'T5: {len(T5)}')
    print(f'Time taken: {time()-start}')
    k0 = []
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    k5 = []
    for op in T5:
        k = op.k
        if k == 0:
            k0.append(op)
        elif k == 1:
            k1.append(op)
        elif k == 2:
            k2.append(op)
        elif k == 3:
            k3.append(op)
        elif k == 4:
            k4.append(op)
        elif k == 5:
            k5.append(op)
    print(f'k0 = {len(k0)}')
    print(f'k1 = {len(k1)}')
    print(f'k2 = {len(k2)}')
    print(f'k3 = {len(k3)}')
    print(f'k4 = {len(k4)}')
    print(f'k5 = {len(k5)}')
    print(f'total = {len(k0)+len(k1)+len(k2)+len(k3)+len(k4)+len(k5)}')
    return T5



T1 = generate_T1(load=True)
T2 = generate_T2(load=True)
T3 = generate_T3()
T4 = generate_T4()
