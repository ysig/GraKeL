cimport cython
cdef extern from "include/functions.hpp":
    unsigned int ArashPartov(const char* str, unsigned int length)
