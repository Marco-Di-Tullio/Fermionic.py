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
        return 0
    else:
        copy[pos][mode] = 0
        new_index = base_list.index(copy[pos]) + 1
        if mode == 0:
            sign = 1
        else:
            sign = (-1)**(sum(copy[pos][0:mode]))
        
        return new_index, sign


# basis, vacio = integer_digits(4)
# assert c(2,basis,2)[0] == 1
# assert c(3,basis,15)[0] == 15
# assert c(1,basis,2) == 0
# assert c(3,basis,3)[1] == -1


def operators(n):
    '''operators generates all the destruction and creation 
    operators in dimension n. It also generates the two order
    operators.
    Ideally, this function would return sparse matrices, but in 
    python I could not find a way of generating a list of sparse
    matrices. In mathematica this can be easily done
    '''
    
    basis, vacio = integer_digits(n)
    l = len(vacio)  
    cm = []
    cd = []
    cdcm = np.zeros((n,n,l,l))
    cmcd = np.zeros((n,n,l,l))
    cmcm = np.zeros((n,n,l,l))
    cdcd = np.zeros((n,n,l,l))
     
    for i in range(n):
        row = []
        col = []
        data = []
        for j in range(l):
            if c(i, basis, j) != 0:
                col.append(j)
                row.append(c(i, basis, j)[0]-1)
                data.append(c(i, basis, j)[1])               
        cmi = sparse.csr_matrix((data, (row, col)), shape=(l, l))
        cdi = sparse.csr_matrix.transpose(cmi)
        cm.append(cmi.toarray())
        cd.append(cdi.toarray())
    
    for i in range(n):
        for j in range(n):
            cdcm[i,j] = cd[i].dot(cm[j])
            cmcd[i,j] = cm[i].dot(cd[j])        
            cdcd[i,j] = cd[i].dot(cd[j])
            cmcm[i,j] = cm[i].dot(cm[j])
    return cm, cd, cdcm, cmcd, cdcd, cmcm

n = int(input('Choose the dimension of the system: '))

cm, cd, cdcm, cmcd, cdcd, cmcm = operators(n)
print("--- %s seconds ---" % (time.time() - start_time))

#def duplicates(lst, item):
#    '''This function is used in c, for finding repeated indexes'''
#    return [i for i, x in enumerate(lst) if x == item]

