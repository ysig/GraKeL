#!/usr/bin/env python
from __future__ import print_function

import sys
from pynauty import autgrp

# List of graphs for testing
#
# Structure:
#   [[name, Graph, numorbit, grpsize, generators]]
#
# numorbit, grpsize, generators was calculated by dreadnut
#

from data_graphs import graphs


if __name__ == '__main__':
    print('Testing pynauty.autgrp()')
    print('Python version: ' + sys.version)
    print('Starting ...')
    passed = 0
    failed = 0
    for gname, g, numorbit, grpsize, gens in graphs:
        print('%-17s ...' % gname, end=' ')
        sys.stdout.flush()
        generators, order, o2, orbits, orbit_no = autgrp(g)
        if generators == gens and orbit_no == numorbit and order == grpsize:
            print('OK')
            passed += 1
        else:
            print('failed')
            failed +=1
    print('... done.')
    if failed > 0:
        print('passed = %d   failed = %d' % (passed, failed))
    else:
        print('All tests passed.')

