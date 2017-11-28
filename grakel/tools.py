""" A python file that implements classes and functions, generally
    useful in other files of the GraKeL project.

"""

from __future__ import generators
import operator

import numpy as np

from scipy.special import binom as binomial

        
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
        value: the value that is going to be added
        keys: a list of keys
    """
    address = dictionary
    for k in keys[:-1]:
        if(k not in address):
            address[k] = dict()
        address = address[k]
    address[keys[-1]] = value
    
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
            if type(k) is list:
                k = tuple(k)
            elif type(k) is set:
                k = frozenset(k)
            else:
                # check for hashability
                # if not raise error
                pass
            if k not in inv:
               inv[k] = list()
            inv[k].append(a)
        # returns a dictionary of lists
    return inv

global ops
ops = {
    '>': operator.gt,
    '<': operator.lt,
    '>=': operator.ge,
    '<=': operator.le,
    '==': operator.eq
    }


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
    if op not in ['>', '<', '>=', '<=', '==']:
        # Raise exception unsupported operator?
        pass

    opr = ops[op]
    
    if(s==-1):
        s = matrix.shape[0]
    
    dictionary = dict()
    for i in range(0,s):
        line = matrix[i,:]
        if not allow_diagonal:
            np.delete(line,i)
        w = np.where(opr(line,const_value))
        dictionary[i] = list(w[0])
    return dictionary

def extract_matrix(mat, a, b):
    """ Subtract the matrix corresponding to two 
        index lists from an original matrix
        
        a, b: the two corresponding index lists
    """
    if (len(a) != len(b)):
        # Raise exception: Index lists must have the same size
        pass
        
    mat_a, mat_b = mat[a,:], mat[b,:]
    
    A = np.concatenate((mat_a[:,a], mat_a[:,b]), axis=1)
    B = np.concatenate((mat_b[:,a], mat_b[:,b]), axis=1)
    
    return np.concatenate((A,B),axis=0)
    
def distribute_samples(n, subsets_size_range, n_samples):
        """ A function that is used in order to 
            distribute evenly, the amount of samples
            that will be drawn from a range of subset
            size, from an original set of given size.
            
            n: set size
            subsets_size_range: A touple having the min and the max subset size
            n_samples: the number of samples
        """
        ## Check input
        min_ss, max_ss = subsets_size_range[0], subsets_size_range[1]
        
        if min_ss <= 1:
            # raise minimum subset size must be bigger than one?
            pass
        if min_ss > max_ss:
            # raise minimum subset size must be smaller than maximum?
            pass
        if min_ss > n:
            # raise minimum subset size must not exceed graph size?
            pass
        if max_ss > n:
            # raise warning maximum subset to big. saving min?
            max_ss = n
            pass
        
        # Distribute samples to subset groups
        availabilities_on_subsets = sorted([(k,int(binomial(n,k))) for k in range(min_ss,max_ss+1)], key = lambda x: x[1])
        n_availabilities = sum(item[1] for item in availabilities_on_subsets)
        
        # Semantic Exception
        if n_availabilities < n_samples:
            # raise error must provide a smaller number of samples?
            pass
        
        samples_on_subsets = dict()
        available_samples = n_samples
        # a variable that helps distributing equally
        cache = 0
        for (i,n) in availabilities_on_subsets:
            a = round((n/n_availabilities)*n_samples)
            value = -1
            if a < n:
                if a > 0:
                    q = a + cache - n
                    if q >= 0:
                        cache = q                    
                        value = n
                    else:
                        value = a + cache
            elif (a >= n):
                cache += a-n
                value = n
                
            # If maximum number of samples is reached break
            if value >= available_samples:
                samples_on_subsets[min_ss+i] = available_samples
                break
            elif value != -1:
                samples_on_subsets[min_ss+i] = value
               
        return samples_on_subsets

