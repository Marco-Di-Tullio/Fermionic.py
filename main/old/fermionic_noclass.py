#fermionic.py

import time
import numpy as np
from scipy import sparse

#start_time = time.time()

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

    rowb = []
    colb = []
    data = []
    vacio = [0 for i in range(2**n)]
    vacio[1] = 1
    for i in range(1,2**n):
        binary_base = format(i, 'b')
        bin_vector = splitter(binary_base)
        while len(bin_vector) < n:
            bin_vector = np.insert(bin_vector,0,0)
        for j in range(n):
            if bin_vector[j] == 1:
                rowb.append(i+1)
                colb.append(j)
                data.append(1)
    basis = sparse.csr_matrix((data, (rowb, colb)))
    return rowb, colb, basis, vacio

def operators(n):
    '''operators generates all the destruction and creation
    operators in dimension n stacked in a sparse matrix
    '''

    rowb, colb, basis, vacio = integer_digits(n)
    l = 2**n
    lb = len(rowb)
    row = []
    col = []
    data = []
    for i in range(lb):
        j = rowb[i]-2**(n-1-colb[i])
        sign = (-1)**(basis[rowb[i],0:colb[i]].sum())
        row.append(j-1)
        col.append(rowb[i]-1+colb[i]*l)
        data.append(sign)
    cm_tot = sparse.csr_matrix((data, (row, col)),shape=(l, n*l))
    cd_tot = sparse.csr_matrix.transpose(cm_tot)
    return cm_tot, cd_tot, l

'''
----------------Operators----------------------

I can write this part in a different file.
When importing an external module, jupyter requests a
different addres than regular python interpreters

If you would like to use this in a separate file:
for running this script in jupyter, use main.fermionic
for runnin this script in regular python, use fermionic

from main.fermionic import *

'''





n = int(input('Choose the dimension of the system: '))
start_time = time.time()
cm_tot, cd_tot, l = operators(n)
print("--- %s seconds ---" % (time.time() - start_time))

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

'''
----------------Class state----------------------
I can write this part in a different file.
When importing an external module, jupyter requests a
different addres than regular python interpreters

If you would like to use this in a separate file:
for running this script in jupyter, use main.fermionic
for runnin this script in regular python, use fermionic

from main.operators import *
'''

class State(object):
    def __init__(self, state):
        self.state = state
        n = int(np.log2(len(self.state)))
        self.dimension = n

    def rhosp(self):
        '''rhosp is the one body density matrix
        '''

        n = self.dimension
        state = self.state
        rhosp = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                rhosp[i,j] = state.dot(cdcm(j,i).dot(state))
        return rhosp

    def eigen(self):
        '''eigen are the eigenvalues of the one body density matrix
        '''
        eigen = np.linalg.eigvalsh(self.rhosp())
        return eigen

    def ssp(self):
        '''ssp is the one body entropy
        '''
        s = -1*np.sum(self.eigen()*np.log(self.eigen())+(1-self.eigen())*np.log(1-self.eigen()))
        return s
