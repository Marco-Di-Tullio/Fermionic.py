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



start_time = time.time()

n = int(input('Choose the dimension of the system: '))

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
            
