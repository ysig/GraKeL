/*
   nautywrap.h

Copyright (c) 2015 Peter Dobsan

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.  This program is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
*/


#include <nauty.h>

#define WORKSPACE_FACTOR    66
#define NUM_GENS_INCR      500

//  a compound data structure to hold all the nauty data structures
//  which describe/used for computing with a given graph

typedef struct {
    optionblk *options;
    
    int         no_vertices;
    int         no_setwords;
    // adjacency matrix as a bit-array
    setword     *matrix;
    // adjacency matrix for the canonical graph
    setword     *cmatrix;
    // coloring: represented as 0-level partition of vertices
    int         *lab;
    int         *ptn;
   // orbits under Autgrp
    int         *orbits;

    // list of generators of Autgrp
    int max_no_generators;
    int no_generators;
    permutation **generator;

    statsblk    *stats;
    int         worksize;
    setword     *workspace;
} NyGraph;

