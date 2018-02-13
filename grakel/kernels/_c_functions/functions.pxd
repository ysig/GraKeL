cimport cython
cdef extern from "include/functions.hpp":
    unsigned int ArashPartov(const char* str, unsigned int length)
    void sm_core_init(double value, int* d, int nv, int kappa, double *cost_vertices, double **cost_edges, double *total_value)
