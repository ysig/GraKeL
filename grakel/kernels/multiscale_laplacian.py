""" The Multiscale Laplacian Graph Kernel as defined
    by [Risi Kondor, Horace Pen (2016)]
"""
import numpy as np

from numpy.linalg import eig, inv, det

from ..graph import graph, laplacian
from ..tools import extract_matrix

def multiscale_laplacian(X, Y, Phi_x, Phi_y, L=3, gamma=0.01):
    """ The Laplacian Graph Kernel
        as proposed by Risi Kondor, 
        Horace Pen at 2016
        
        X,Y: Valid graph formats to be compared
        Phi_{x,y}: The corresponding feature values dictionaries
                    
        gamma: A small softening parameter of float value
    """
    Gx = graph(X, Phi_x)
    Gy = graph(Y, Phi_y)
    return multiscale_laplacian_inner(Gx, Gy, gamma=gamma)
    
    
def multiscale_laplacian_inner(Gx, Gy, L=3, gamma=0.01):
    """ The Laplacian Graph Kernel
        as proposed by Risi Kondor, 
        Horace Pen at 2016
        
        Gx,Gy: Valid graph formats having as labels the phi_x, phi_y
        
        gamma: A small softening parameter of float value
        L: number of neighbourhoods
    """
    
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
    
    # calculate the neighboorhoods
    Nx = make_nested_neighbourhoods(Gx, L=L)
    Ny = make_nested_neighbourhoods(Gy, L=L)
    
    # a lambda that calculates indexes inside the gram matrix
    # and the corresponindg laplacian given a node and a level
    pick = lambda node, level: (Nx[node][level], laplacian(Ax[Nx[node][level],:][:,Nx[node][level]])) if node<Gx.n else ([idx+Gx.n for idx in Ny[node-Gx.n][level]], laplacian(Ay[Ny[node-Gx.n][level],:][:,Ny[node-Gx.n][level]]))

    for l in range(1,L+1):
        gm = gram_matrix
        gram_matrix = np.empty(shape=(gram_matrix_size,gram_matrix_size))
        for i in range(0, gram_matrix_size):
            # calculate the correct indexes of neighbours
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
    """ Helping function for multiscale gaussian
       L_{x,y}: Laplacians of graph {x,y}
                numpy arrays of size n_{x,y} times n_{x,y}
       Gram_matrix: The corresponding gram matrix for the two graphs
                    numpy array of size n_x + n_y
       gamma: A smoothness parameter [float] - close to zero.
       
    """
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

def make_nested_neighbourhoods(G, L):
    """ Calculates nested neighborhoods needed
        for multiscale laplacian calculation
        
        G: a graph type object
        L: number of neighbourhoods for nesting
    """
    N = dict()
    
    # initialization
    for i in range(0,G.n):
        N[i] = dict()
        N[i][1] = sorted([i]+list(G.neighbours(i)))
        
    # calculate neighbourhoods
    # by a recursive formula
    # for all levels from 2 to L
    for level in range(1,L):
        for i in range(0,G.n):
            neighbours = set()
            for w in N[i][level]:
                neighbours = neighbours.union(set(N[w][level]))
            N[i][level+1] = sorted(list(neighbours))
    
    return N
