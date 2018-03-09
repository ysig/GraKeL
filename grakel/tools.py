"""A file having general functions and classes usefull insid `grakel`."""

from __future__ import generators
import collections
import operator
import warnings

import numpy as np

from scipy.special import binom as binomial


class priority_dict(dict, object):
    """A priority dictionary.

    This implementation of priority dictionaries using binary heaps
    was implemented by David Eppstein, UC Irvine, 8 Mar 2002
    and was found at:
    http://code.activestate.com/recipes/117228-priority-dictionary/
    """

    def __init__(self):
        """Initialize a priority heap object."""
        self.__heap = []
        dict.__init__(self)

    def smallest(self):
        """Find smallest item after removing deleted items from heap."""
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
        """Create destructive sorted iterator of priorityDictionary."""
        def iterfn():
            while len(self) > 0:
                x = self.smallest()
                yield x
                del self[x]
        return iterfn()

    def __setitem__(self, key, val):
        """`__setitem__` default method for primary heap.

        Change value stored in dictionary and add corresponding
        pair to heap.  Rebuilds the heap if the number of deleted items grows
        too large, to avoid memory leakage.
        """
        dict.__setitem__(self, key, val)
        heap = self.__heap
        if len(heap) > 2 * len(self):
            self.__heap = [(v, k) for k, v in self.iteritems()]
            self.__heap.sort()  # builtin sort likely faster than O(n) heapify
        else:
            newPair = (val, key)
            insertionPoint = len(heap)
            heap.append(None)
            while insertionPoint > 0 and \
                    newPair < heap[(insertionPoint-1)//2]:
                heap[insertionPoint] = heap[(insertionPoint-1)//2]
                insertionPoint = (insertionPoint-1)//2
            heap[insertionPoint] = newPair

    def setdefault(self, key, val):
        """Reimplement setdefault to call our customized __setitem__."""
        if key not in self:
            self[key] = val
        return self[key]


def nested_dict_add(dictionary, value, *keys):
    """Nested adding in a dictionary.

    Parameters
    ----------
        dictionary : dict
            The dictionary that is going to be parsed.

        value : object
            The value that is going to be added.

        keys : list(hashable)
            A list of keys.
    Returns
    -------
    None.

    """
    address = dictionary
    for k in keys[:-1]:
        if(k not in address):
            address[k] = dict()
        address = address[k]
    address[keys[-1]] = value


def nested_dict_get(dictionary, *keys, **kargs):
    """Get an item from a nested dictionary.

    Checks if an item exists in a nested level by keys in a dictionary and if
    yes returns it. Otherwise return default.

    Parameters
    ----------
        dictionary : dict
            The dictionary that is going to be parsed.

        keys : list(hashable)
            A list of keys.

        default : object, default=None
            The value to return, if the element is not found.

    Returns
    -------
        dict_elem : object
            Returns either the dictionary element or the value in default.

    """
    if len(kargs) == 1 and "default":
        default = kargs["default"]
    elif len(kargs) == 0:
        default = None
    else:
        raise TypeError('optional argument can only be "default"')

    element = dictionary
    for k in keys:
        if (k in element):
            element = element[k]
        else:
            return default
    return element


def inv_dict(d):
    """Calculate the inverse dictionary.

    Parameters
    ----------
        d: dict
            A `surjective` dictionary.

    Returns
    -------
        inv_dict : dict
            The inverse dictionary.

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
                if not isinstance(k, collections.Hashable):
                    raise ValueError('in order to calculate inverse \
                          dictionary, values must be hashable')
            if k not in inv:
                inv[k] = list()
            inv[k].append(a)
        # returns a dictionary of lists
    return inv


ops = {
    '>': operator.gt,
    '<': operator.lt,
    '>=': operator.ge,
    '<=': operator.le,
    '==': operator.eq
    }


def matrix_to_dict(matrix, op='==', const_value=0, allow_diagonal=False):
    """Matrix to dictionary transformation based on constraints.

    Parameters
    ----------
        matrix : np.array
            A square to be transformed matrix.
        op : str
            An operator applied as a constraint.
        const_value : object
            The value to be compared with.
        allow_diagonal: bool
            Allows matrix diagonal to be added at input.

    Returns
    -------
        produced_dict : dict
            The produced matrix dictionary.

    """
    if op not in ['>', '<', '>=', '<=', '==']:
        raise ValueError('unsupported operator')

    opr = ops[op]

    s = matrix.shape[0]
    dictionary = dict()
    for i in range(0, s):
        line = matrix[i, :]
        if not allow_diagonal:
            np.delete(line, i)
        w = np.where(opr(line, const_value))
        dictionary[i] = list(w[0])
    return dictionary


def distribute_samples(n, subsets_size_range, n_samples):
    """Distribute samples evenly in a given range.

    A function that is used in order to distribute evenly, the amount of
    samples that will be drawn from a range of subset's sizes, from an
    original set of given size.

    Parameters
    ----------
        n : int
            The set size.

        subsets_size_range : tuple
            A touple having the min and the max subset size.

        n_samples : int
            The number of samples.

    Returns
    ------
        samples_on_subsets : dict
            Returns a dictionary of samples, for each subset.

    """
    # Check input
    min_ss, max_ss = subsets_size_range[0], subsets_size_range[1]

    if min_ss <= 1:
        raise ValueError('minimum subset size must be bigger than one')
    if min_ss > max_ss:
        raise ValueError('minimum subset size must be \
        smaller than maximum')
    if min_ss > n:
        raise ValueError('minimum subset size must not exceed graph size')
    if max_ss > n:
        warnings.warn('maximum subset size to big - adjusting to set size')
        max_ss = n

    # Distribute samples to subset groups
    availabilities_on_subsets = sorted([(k, int(binomial(n, k)))
                                        for k in range(min_ss, max_ss+1)],
                                       key=lambda x: x[1])

    n_availabilities = sum(item[1] for item in availabilities_on_subsets)

    # Semantic Exception
    if n_availabilities < n_samples:
        warnings.warn(
            'number of samples exceedes the number of availabilities, ' +
            'n_samples is now replaced by the availabilities_on_subsets')
        return dict(availabilities_on_subsets)

    samples_on_subsets = dict()
    available_samples = n_samples
    # a variable that helps distributing equally
    cache = 0
    for (i, n) in availabilities_on_subsets:
        a = int(round((n/(n_availabilities*1.0))*n_samples))
        value = -1
        if a < n:
            if a > 0:
                q = a + cache - n
                if q >= 0:
                    cache = q
                    value = n
                else:
                    value = a + cache
        elif a >= n:
            cache += a-n
            value = n

        # If maximum number of samples is reached break
        if value >= available_samples:
            samples_on_subsets[min_ss+i] = available_samples
            break
        elif value != -1:
            samples_on_subsets[min_ss+i] = value

    return samples_on_subsets
