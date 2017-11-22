""" A python file that implements classes and functions, generally
    useful in other files of the GraKeL project.

"""

from __future__ import generators
import operator

import numpy as np
        
class priority_dict(dict,object):
    def __init__(self):
        """Initialize priorityDictionary by creating binary heap of pairs (value,key).
           Note that changing or removing a dict entry will
           not remove the old pair from the heap until it is found by smallest() or
           until the heap is rebuilt.

           This implementation of priority dictionaries using binary heaps
           was implemented by David Eppstein, UC Irvine, 8 Mar 2002
           and was found at: http://code.activestate.com/recipes/117228-priority-dictionary/
        """
        self.__heap = []
        dict.__init__(self)


    def smallest(self):
        """ Find smallest item after removing deleted items from heap. """
        if len(self) == 0:
            raise(IndexError, "smallest of empty priorityDictionary")
        heap = self.__heap
        while heap[0][1] not in self or self[heap[0][1]] != heap[0][0]:
            lastItem = heap.pop()
            insertionPoint = 0
            while 1:
                smallChild = 2*insertionPoint+1
                if smallChild+1 < len(heap) and \
                        heap[smallChild] > heap[smallChild+1]:
                    smallChild += 1
                if smallChild >= len(heap) or lastItem <= heap[smallChild]:
                    heap[insertionPoint] = lastItem
                    break
                heap[insertionPoint] = heap[smallChild]
                insertionPoint = smallChild
        return heap[0][1]
    
    def __iter__(self):
        """ Create destructive sorted iterator of priorityDictionary.
        
        """
        def iterfn():
            while len(self) > 0:
                x = self.smallest()
                yield x
                del self[x]
        return iterfn()
    
    def __setitem__(self,key,val):
        """ Change value stored in dictionary and add corresponding
            pair to heap.  Rebuilds the heap if the number of deleted items grows
            too large, to avoid memory leakage.

        """
        dict.__setitem__(self,key,val)
        heap = self.__heap
        if len(heap) > 2 * len(self):
            self.__heap = [(v,k) for k,v in self.iteritems()]
            self.__heap.sort()  # builtin sort likely faster than O(n) heapify
        else:
            newPair = (val,key)
            insertionPoint = len(heap)
            heap.append(None)
            while insertionPoint > 0 and \
                    newPair < heap[(insertionPoint-1)//2]:
                heap[insertionPoint] = heap[(insertionPoint-1)//2]
                insertionPoint = (insertionPoint-1)//2
            heap[insertionPoint] = newPair
    
    def setdefault(self,key,val):
        """ Reimplement setdefault to call our customized __setitem__.

        """
        if key not in self:
            self[key] = val
        return self[key]

def nested_dict_add(dictionary, value, *keys):
    """ Does a nesting adding in a dictionary
        dictionary: the dictionary that is going to be parsed

        dictionary: the dictionary that is going to be parsed
        value: the value that is going to be added
        keys: a list of keys
    """
    address = dictionary
    l = len(keys)
    i = 1
    for k in keys:
        if(i<l):
            if(k not in address):
                address[k] = dict()
            address = address[k]
        elif(i==l):
            address[k] = value            
        i+=1
        
def nested_dict_get(dictionary, *keys):
    """ Checks if an item exists in a nested
        level by keys in a dictionary and if
        yes returns it. Otherwise return None

        dictionary: the dictionary that is going to be parsed
        keys: a list of keys
        """
    element = dictionary
    for k in keys:
        if (k in element):
            element = element[k]
        else:
            return None
    return element
    

def inv_dict(d):
    """ A function that given a surjective
        dictionary, calculates its dictionary inverse
            
        d: dictionary
    """
    inv = dict()
    if bool(d):
        for a in d.keys():
            k = d[a]
            if k not in inv:
               inv[k] = list()
            inv[k].append(a)
        # returns a dictionary of list items
    return inv

def matrix_to_dict(matrix, op='==', const_value=0, s=-1, allow_diagonal=False):
    """ Transform a matrix to a dictionary
        adding for each row only the elements
        that satisfy as certain constraint
        specifically.

        matrix: A square numpy matrix
        op: an operator applied as a constraint
        const_value: the value to be compared with
        s: if -1 calculates the matrix shape
           else takes it as input
        allow_diagonal: allows matrix diagonal to be added
                        at input
    """
    ops = {
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        '==': operator.eq
        }
    opr = ops[op]
    if(s==-1):
        s = matrix.shape[0]
    dictionary = dict()
    for i in range(0,s):
        for j in range(0,s):
            if(allow_diagonal or i!=j):
                if(opr(matrix[i,j],const_value)):
                    if i not in dictionary:
                        dictionary[i] = list()
                    dictionary[i].append(j)
    return dictionary

def extract_matrix(mat, a, b):
    """ Subtract the matrix corresponding to two 
        index lists from an original matrix
        
        a, b: the two corresponding index lists
    """
    mat_a, mat_b = mat[a,:], mat[b,:]
    
    A = np.concatenate((mat_a[:,a], mat_a[:,b]), axis=1)
    B = np.concatenate((mat_b[:,a], mat_b[:,b]), axis=1)
    
    return np.concatenate((A,B),axis=0)
