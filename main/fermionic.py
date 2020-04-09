#fermionic.py

import numpy as np
from scipy import sparse
import time


start_time = time.time()

def splitter(x):
    '''This function creates an numpy array where each 
    element is a digit from the original number
    '''
    
    vectorized = []
    for i in str(x):
        vectorized.append(int(i))
    vectorized = np.array(vectorized)
    return vectorized
        

def integer_digits(n):
    '''integer_digits(n) creates a fermionic basis for dimension n
    '''
    
    assert n > 0, 'n must be greater than 0'
    basis = np.array([0 for i in range(n)])
    for i in range(1,2**n):
        binary_base = format(i, 'b')
        bin_vector = splitter(binary_base)
        while len(bin_vector) < n:
            bin_vector = np.insert(bin_vector,0,0)
        basis = np.vstack((basis,bin_vector))
        vacio = np.zeros(len(basis))
        vacio[0] = 1
    return basis, vacio

def c(mode, base, pos):
    '''c(mode, base, pos) gives the index in the base
    of the element pos after the application of the 
    destruction operator over mode
    '''
    
    base_list = base.tolist()
    copy = base_list
    if copy[pos][mode] == 0:
        return None
    else:
        copy[pos][mode] = 0
        new_index = base_list.index(copy[pos]) + 1
        if mode == 0:
            sign = 1
        else:
            sign = (-1)**(sum(copy[pos][0:mode]))
        
        return new_index, sign

def operators(n):
    '''operators generates all the destruction and creation 
    operators in dimension n stacked in a sparse matrix
    '''
    
    basis, vacio = integer_digits(n)
    l = len(vacio)   
    row = []
    col = []
    data = []    
    for i in range(n):
        for j in range(l):
            if c(i, basis, j):
                row.append(c(i, basis, j)[0]-1)
                col.append(j+i*l)
                data.append(c(i, basis, j)[1])               
    cm_tot = sparse.csr_matrix((data, (row, col)), shape=(l, n*l))
    cd_tot = sparse.csr_matrix.transpose(cm_tot)              
    return cm_tot, cd_tot, l


n = int(input('Choose the dimension of the system: '))

cm_tot, cd_tot, l = operators(n)

def cm(i, cm_tot = cm_tot, l = l):
    '''cm(i) cuts the big cm_tot sparse matrix to the
    corresponding mode i
    '''
    
    return cm_tot[:,l*i:l*(i+1)]

def cd(i, cd_tot = cd_tot, l = l):
    '''cd(i) cuts the big cd_tot sparse matrix to the
    corresponding mode i
    '''
    
    return cd_tot[l*i:l*(i+1),:]
    
def cdcm(i,j):    
    cdi = cd(i)
    cmj = cm(j)
    return cdi.dot(cmj)

def cmcd(i,j):    
    cmi = cm(i)
    cdj = cd(j)
    return cmi.dot(cdj)

def cmcm(i,j):    
    cmi = cm(i)
    cmj = cm(j)
    return cmi.dot(cmj)

def cdcd(i,j):    
    cdi = cd(i)
    cdj = cd(j)
    return cdi.dot(cdj)



print("--- %s seconds ---" % (time.time() - start_time))

#def duplicates(lst, item):
#    '''This function is used in c, for finding repeated indexes'''
#    return [i for i, x in enumerate(lst) if x == item]


# def operators(n):
#     '''operators generates all the destruction and creation 
#     operators in dimension n. It also generates the two order
#     operators.
#     Ideally, this function would return sparse matrices, but in 
#     python I could not find a way of generating a list of sparse
#     matrices. In mathematica this can be easily done
#     '''
    
#     basis, vacio = integer_digits(n)
#     l = len(vacio)  
#     cm = []
#     cd = []

     
#     for i in range(n):
#         row = []
#         col = []
#         data = []
#         for j in range(l):
#             if c(i, basis, j) != 0:
#                 col.append(j)
#                 row.append(c(i, basis, j)[0]-1)
#                 data.append(c(i, basis, j)[1])               
#         cmi = sparse.csr_matrix((data, (row, col)), shape=(l, l))
#         cdi = sparse.csr_matrix.transpose(cmi)
#         cm.append(cmi.toarray())
#         cd.append(cdi.toarray())
    
    
#     return cm, cd

# def two_body_op(n, cm_tot = cm_tot, cd_tot = cd_tot):
#     '''This function defines the two body operators
#     '''
#     basis, vacio = integer_digits(n)
#     l = len(vacio)  
#     cdcm = np.zeros((n,n,l,l))
#     cmcd = np.zeros((n,n,l,l))
#     cmcm = np.zeros((n,n,l,l))
#     cdcd = np.zeros((n,n,l,l))
#     for i in range(n):
#         for j in range(n):
#             cdcm[i,j] = cd_tot[i].dot(cm_tot[j])
#             cmcd[i,j] = cm_tot[i].dot(cd_tot[j])        
#             cdcd[i,j] = cd_tot[i].dot(cd_tot[j])
#             cmcm[i,j] = cm_tot[i].dot(cm_tot[j])
#     return cdcm, cmcd, cdcd, cmcm

# from sympy.physics.quantum import AntiCommutator, Commutator
# from sympy.physics.quantum.fermion import FermionOp
# from sympy.physics.quantum import Operator

# def com(A,B):
#     A = Operator('A')
#     B = Operator('B')
#     return Commutator(A,B).expand(commutator=True)