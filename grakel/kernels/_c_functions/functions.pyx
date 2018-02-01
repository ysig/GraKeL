from libc.string cimport const_char
cimport functions

def APHash(word):
    """C++ wrapped implementation of Arash Partov Hashing."""
    bs = word.encode('UTF-8')
    cdef int length = len(word);
    cdef const_char* string = bs;
    return functions.ArashPartov(string, length)
