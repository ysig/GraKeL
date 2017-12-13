""" The Multiscale Laplacian Graph Kernel as defined in :cite:`Kondor2016TheML`
"""
import warnings

import numpy as np

from numpy.linalg import eig, inv, det

from ..graph import graph, laplacian
from ..tools import extract_matrix

def multiscale_laplacian(X, Y, Phi_x, Phi_y, L=3, gamma=0.01):
    """ The Laplacian Graph Kernel as proposed in :cite:`Kondor2016TheML`
        
    arguments:    
        - X,Y (valid graph format): the pair of graphs on which the kernel is applied
        - Phi_{x,y} (dict): corresponding feature vectors for nodes
        - gamma (float): A small softening parameter of float value
        - L (int): number of neighborhoods
        
    returns:
        number. The kernel value        
    """
    Gx = graph(X, Phi_x)
    Gy = graph(Y, Phi_y)
    return multiscale_laplacian_inner(Gx, Gy, gamma=gamma)
    
    
def multiscale_laplacian_inner(Gx, Gy, L=3, gamma=0.01):
    """ The Laplacian Graph Kernel as proposed in :cite:`Kondor2016TheML`
    
    arguments:
        - Gx, Gy (graph): the pair of graphs on which the kernel is applied with feature vectores as node labels
        - gamma (float): a small softening parameter
        - L (int): number of neighborhoods

    returns:
        number. The kernel value
    """
    if gamma > 0.05:
        warnings.warn('gamma to big')
    elif gamma == .0:
        warnings.warn('with zero gamma the calculation may crash')
    elif gamma < 0:
        raise ValueError('gamma must be positive')
    
    if type(L) is not int:
        raise ValueError('L must be an integer')
    elif L < 0:
        raise ValueError('L must be positive')
    
    # Set desired format
    Gx.desired_format("adjacency")
    Gy.desired_format("adjacency")
    
    # calculate adjacency matrix
    Ax = Gx.adjacency_matrix
    Ay = Gy.adjacency_matrix
    
    # Get Feature Vectors
    phi_x = Gx.get_labels()
    phi_y = Gy.get_labels()
    
    # create the gram matrix
    gram_matrix_size = Gx.n + Gy.n
    gram_matrix = np.empty(shape=(gram_matrix_size,gram_matrix_size))
    
    # initialise the gram matrix
    pick = lambda i: list(phi_x[i]) if i<Gx.n else list(phi_y[i-Gx.n])
    for i in range(0, gram_matrix_size):
        vec_a = pick(i)
        for j in range(0, gram_matrix_size):
            vec_b = pick(i)
            gram_matrix[i,j] = np.dot(vec_a, vec_b)
    
    # calculate the neighborhoods
    Nx = Gx.produce_neighborhoods(r=L)
    Ny = Gy.produce_neighborhoods(r=L)
    
    # a lambda that calculates indexes inside the gram matrix
    # and the corresponindg laplacian given a node and a level
    pick = lambda node, level: (Nx[level][node], laplacian(Ax[Nx[level][node],:][:,Nx[level][node]])) if node<Gx.n else ([idx+Gx.n for idx in Ny[level][node-Gx.n]], laplacian(Ay[Ny[level][node-Gx.n],:][:,Ny[level][node-Gx.n]]))

    for l in range(1,L+1):
        gm = gram_matrix
        gram_matrix = np.empty(shape=(gram_matrix_size,gram_matrix_size))
        for i in range(0, gram_matrix_size):
            # calculate the correct indexes of neighbors
            # and the corresponding laplacian
            (idx_i, La) = pick(i,l)
            for j in range(0, gram_matrix_size):
                (idx_j, Lb) = pick(j,l)
                    
                # calculate the corresponding gram matrix
                extracted_gm = extract_matrix(gram_matrix, idx_i, idx_j)
                
                # calculate the core for this step
                gram_matrix[i,j] = generalized_FLG_core(La, Lb, extracted_gm, gamma=gamma)
    
    # Calculate the full laplacian
    Lx = Gx.laplacian()
    Ly = Gy.laplacian()
            
    return generalized_FLG_core(Lx, Ly, gram_matrix, gamma=gamma)
    
def generalized_FLG_core(Lx, Ly, gram_matrix, gamma=0.05, heta=0.1):
    """ Helping function for the multiscale gaussian
    
    arguments:
        - L_{x,y} (np.array): Laplacians of graph {x,y}
        - Gram_matrix (np.array): The corresponding gram matrix for the two graphs
        - gamma (float): a small softening parameter
        - heta (float): number of neighborhoods

    returns:
        number. The FLG core kernel value
    """
    if heta > 0.2:
        warnings.warn('heta to big')
    elif gamma == .0:
        warnings.warn('with zero heta the calculation may crash')
    elif heta < 0:
        raise ValueError('heta must be positive')
 
 
    nx = Lx.shape[0]
    ny = Ly.shape[0]
    
    # w the eigen vectors, v the eigenvalues
    v, w = eig(gram_matrix)
    n_eig = len(v)

    # keep only the positive
    vp = {i: v[i] for i in range(0,len(v)) if (v[i]>.0)}
    valid_keys = sorted(list(vp.keys()))
    n_peig = len(valid_keys)

    # calculate the Q matrix
    Q = np.empty(shape=(n_eig,n_peig))
    for (k,idx) in zip(vp.keys(),range(0,n_peig)):
        Q[:,idx] = np.sqrt(vp[k]) * w[:,k]

    Qx = Q[0:nx,:]
    Qy = Q[nx:,:]
    
    # Calculate the S matrices
    Tx = np.dot(np.dot(Qx.T,inv(np.add(Lx,heta*np.eye(Lx.shape[0])))),Qx)
    Ty = np.dot(np.dot(Qy.T,inv(np.add(Ly,heta*np.eye(Ly.shape[0])))),Qy)

    Sx = Tx + gamma*np.eye(Tx.shape[0])
    Sy = Ty + gamma*np.eye(Ty.shape[0])
    
    # A small lambda to calculate ^ 1/4
    quatre = lambda x: np.sqrt(np.sqrt(x))
    
    ### !!! Overflow problem!: det(inv(Sx)+inv(Sy))=0.0
    ### Need to solve. Not all executions are fatal
    # Caclulate the kernel nominator
    k_nom = np.sqrt(det(abs(np.multiply(inv(inv(Sx)+inv(Sy)),2))))
    
    # Caclulate the kernel denominator
    k_denom = quatre(abs(det(Sx)*det(Sy)))
    
    return k_nom/k_denom
