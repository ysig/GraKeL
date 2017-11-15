/*
   nautywrap.c

Copyright (c) 2015 Peter Dobsan

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.  This program is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
*/

#include <Python.h>
#include <nauty.h>
#include <nautywrap.h>


//  static global (yuck) variables  -------------------------------------------

// Nauty default options
// unfortunately, defined as macro in Nauty
//
static DEFAULTOPTIONS(default_options);

// a static global handle the current NyGraph object
// needed since there is no way to pass this parameter to store_generator()
static NyGraph *GRAPH_PTR;

//  Utilities  ================================================================

static void store_generator(int count,
        permutation *perm,
        int *orbits,
        int no_orbits,
        int stabvertex,
        int n)
// this function is called by nauty every time a new generator
// of the automorphismgroup of the graph found.
{
    NyGraph *g = GRAPH_PTR;
    int i;
    permutation *p;
    permutation **new;

    if ((p = malloc(n * sizeof(permutation))) == NULL) {
        fprintf(stderr, "Failed to allocate memory for generator #%d.\n",
                g->no_generators);
        exit(1);
    }

    for (i=0; i < n; i++) {
        p[i] = perm[i];
    }
    g->generator[g->no_generators] = p;
    g->no_generators++;

    if (g->no_generators >= g->max_no_generators) {
        // allocate a larger array
        g->max_no_generators += NUM_GENS_INCR;
        new=(permutation **)malloc(g->max_no_generators * sizeof(permutation*));
        if (new == NULL) {
            fprintf(stderr, "Failed to allocate extension for generators.\n");
            exit(1);
        }
        // copy over the old one
        for (i = 0; i < g->no_generators; i++) {
            new[i] = g->generator[i];
        }
        // free the old generator array
        free(g->generator);
        // switch to the new one
        g->generator = new;
    }
}


void destroy_nygraph(NyGraph *g)
//
// free all the allocated memeory for NyGraph
//
{
    int i;

    free(g->options);
    free(g->matrix);
    free(g->cmatrix);
    free(g->lab);
    free(g->ptn);
    free(g->orbits);
    free(g->stats);
    free(g->workspace);

    for (i = 0; i < g->no_generators; i++) {
        free(g->generator[i]);
    }

    free(g);
}


NyGraph * create_nygraph(int no_vertices)
//
// Allocate a data structure to hold all data structures used by Nauty
// 
{
    NyGraph *g;
    int i, *p;

    if (no_vertices < 0) return NULL;

    if ((g = malloc(sizeof(NyGraph))) == NULL) return NULL;

    // initialize evrything here so if allocation fails at any stage later
    // we can clean up
    g->matrix = g->cmatrix = NULL;
    g->lab = g->ptn = g->orbits = NULL;
    g->options = NULL;
    g->stats = NULL;
    g->workspace = NULL;
    g->max_no_generators = 0;
    g->no_generators = 0;
    g->generator = NULL;

    g->no_vertices = no_vertices;
    g->no_setwords = (no_vertices + WORDSIZE - 1) / WORDSIZE;
    nauty_check(WORDSIZE, g->no_setwords, g->no_vertices, NAUTYVERSIONID);

    if ((g->matrix = malloc((size_t) g->no_setwords * 
                    (size_t) (no_vertices * sizeof(setword)))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }
    for (i = 0; i < no_vertices; i++) {
        EMPTYSET((GRAPHROW(g->matrix, i, g->no_setwords)), g->no_setwords);
    }

    g->cmatrix = NULL;

    if ((g->lab = malloc(no_vertices * sizeof(int))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }
    for (i = 0, p = g->lab; i < no_vertices; i++) p[i] = 0;

    if ((g->ptn = malloc(no_vertices * sizeof(int))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }
    for (i = 0, p = g->ptn; i < no_vertices; i++) p[i] = 0;
    
    if ((g->orbits = malloc(no_vertices * sizeof(int))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }
    for (i = 0, p = g->orbits; i < no_vertices; i++) p[i] = 0;

    if ((g->options = malloc(sizeof(optionblk))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }
    memcpy(g->options, &default_options, sizeof(optionblk));

    g->options->digraph = FALSE;
    g->options->getcanon = FALSE;
    g->options->defaultptn = TRUE;      // default is no coloring
    g->options->writeautoms = FALSE;
    g->options->cartesian = TRUE;
    g->options->linelength = 0;
    g->options->userautomproc = store_generator;

    if ((g->stats = malloc(sizeof(statsblk))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }

    g->worksize = WORKSPACE_FACTOR * g->no_setwords;
    if ((g->workspace = malloc(g->worksize * sizeof(setword))) == NULL) {
        destroy_nygraph(g);
        return NULL;
    }

    GRAPH_PTR = g;
    return g;
}


void make_edge(NyGraph *g, int i, int j)
    // connect vertex i to vertex j with an edge
    // if the graph is directed that means the edge (i)----->(j)
{
    ADDELEMENT((GRAPHROW(g->matrix, i, g->no_setwords)), j);
    if (g->options->digraph == FALSE) {
        ADDELEMENT((GRAPHROW(g->matrix, j, g->no_setwords)), i);
    }
}


NyGraph * extend_canonical(NyGraph *g)
// extend the NyGraph structure with cmatrix to hold
// the canonically labeled graph
{
    if ((g->cmatrix = malloc(g->no_setwords * g->no_vertices * WORDSIZE))
            == NULL) {
        destroy_nygraph(g);
        return NULL;
    }
    return g;
}


//  Python functions  =========================================================

static int set_partition(PyObject *py_graph, int *lab, int *ptn)
// Convert the vertex_coloring attribute of a NyGraph object
// into nauty (lab, ptn) data structure at partition level 0
{
    PyObject *partition;
    PyObject *pyset;
    PyObject *iterator;
    PyObject *item;
    int no_parts;
    int i;
    int n;
    int x;

    if ( !(partition = PyObject_GetAttrString(py_graph, "vertex_coloring")) ) {
        PyErr_SetString(PyExc_TypeError,
                "missing 'vertex_coloring' attribute");
        return 0;       // error
    }

    if ((no_parts = PyObject_Length(partition)) <= 0) {
        return -1;      // no coloring
    }

    for (i = n = 0; i < no_parts; i++) {
        pyset = PyList_GET_ITEM(partition, i);
        iterator = PyObject_GetIter(pyset);
        
        while ((item = PyIter_Next(iterator))) {
#if PY_MAJOR_VERSION >= 3
            x = PyLong_AS_LONG(item);
#else
            x = PyInt_AS_LONG(item);
#endif
            Py_DECREF(item);
            lab[n] = x;
            ptn[n++] = 1;
        }
        if (n > 0) {
            ptn[n-1] = 0;
        }

        Py_DECREF(iterator);
    }

    Py_DECREF(partition);

    return 1;           // coloring
}


static PyObject* py_auto_group(NyGraph *g)
// convert generators, orbits etc. into Python representation
// and return it in a tuple:
//      (generators, order, orbits, orbit_no)
{
    int i, j;
    PyObject *py_autgrp;
    PyObject *py_gens;
    PyObject *py_perm;
    PyObject *py_orbits;
    PyObject *py_grpsize1;
    PyObject *py_grpsize2;

    // generators
    py_gens = PyList_New(g->no_generators);
    for (i=0; i < g->no_generators; i++) {
        py_perm = PyList_New(g->no_vertices);
        for (j=0; j < g->no_vertices; j++) {
            PyList_SET_ITEM(py_perm, j,
                    Py_BuildValue("i", g->generator[i][j]));
        }
        PyList_SET_ITEM(py_gens, i, py_perm);
    }

    // group order
    //
    // "the order of the autgrp is equal to
    //      grpsize1 * 10^grpsize2
    // within rounding error" -- Nauty manual
    //
    py_grpsize1 = Py_BuildValue("d", g->stats->grpsize1);
    py_grpsize2 = Py_BuildValue("i", g->stats->grpsize2);

    // orbits
    py_orbits = PyList_New(g->no_vertices);
    for (i=0; i < g->no_vertices; i++) {
        PyList_SET_ITEM(py_orbits, i, Py_BuildValue("i", g->orbits[i]));
    }

    // create return value tuple:
    //      (generators, grpsize1, grpsize2, orbits, orbit_no)
    py_autgrp = PyTuple_New(5);
    PyTuple_SET_ITEM(py_autgrp, 0, py_gens);
    PyTuple_SET_ITEM(py_autgrp, 1, py_grpsize1);
    PyTuple_SET_ITEM(py_autgrp, 2, py_grpsize2);
    PyTuple_SET_ITEM(py_autgrp, 3, py_orbits);
    PyTuple_SET_ITEM(py_autgrp, 4, Py_BuildValue("i", g->stats->numorbits));

    return py_autgrp;
}


NyGraph * _make_nygraph(PyObject *py_graph)
// Convert the Python NyGraph object into a Nauty/C NyGraph object
// and set Nauty options.
{
    NyGraph *g;
    int n;
    set *rowp;
 
    PyObject *adjdict;
    PyObject *key;
    PyObject *adjlist;
    PyObject *p;

    int i;
    int adjlist_length;
    int x, y;

    // get the number of vertices
    if ((p = PyObject_GetAttrString(py_graph, "number_of_vertices")) == NULL) {
        PyErr_SetString(PyExc_TypeError,
                "Missing 'number_of_vertices' attribute");
        return NULL;
    }
#if PY_MAJOR_VERSION >= 3
    n = PyLong_AS_LONG(p);
#else
    n = PyInt_AS_LONG(p);
#endif
    Py_DECREF(p);

    // create an empty Nauty NyGraph object
    if ((g = create_nygraph(n)) == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Nauty NyGraph creation failed");
        return NULL;
    }

    // get directed attribute
    if ((p = PyObject_GetAttrString(py_graph, "directed")) == NULL) {
        PyErr_SetString(PyExc_TypeError, "missing 'directed' attribute");
        return NULL;
    }
    Py_DECREF(p);
    if (PyObject_IsTrue(p)) {
        g->options->digraph = TRUE;
    } else {
        g->options->digraph = FALSE;
    }

    // get the adjacency list dictionary object
    if ((adjdict = PyObject_GetAttrString(py_graph, "adjacency_dict")) == NULL) {
        PyErr_SetString(PyExc_TypeError, "missing 'adjacency_dict' attribute");
        return NULL;
    }

    // iterate over the adjacency list setting
    // the adjacency matrix in the Nauty NyGraph g
    Py_ssize_t pos = 0;
    while (PyDict_Next(adjdict, &pos, &key, &adjlist)) {
#if PY_MAJOR_VERSION >= 3
        x = PyLong_AS_LONG(key);
#else
        x = PyInt_AS_LONG(key);
#endif
        adjlist_length =  PyObject_Length(adjlist);
        rowp = GRAPHROW(g->matrix, x, g->no_setwords);
        for (i=0; i < adjlist_length; i++) {
            p = PyList_GET_ITEM(adjlist, i);
#if PY_MAJOR_VERSION >= 3
            y = PyLong_AS_LONG(p);
#else
            y = PyInt_AS_LONG(p);
#endif
            ADDELEMENT(rowp, y);
            if (g->options->digraph == FALSE) {
                ADDELEMENT((GRAPHROW(g->matrix, y, g->no_setwords)), x);
            }
        }
    }

    Py_DECREF(adjdict);

    // take care of coloring
    x = set_partition(py_graph, g->lab, g->ptn);
    if (x < 0) {
        g->options->defaultptn = TRUE;
    } else if (x == 0) {
        return 0;
    } else {
        g->options->defaultptn = FALSE;
    }

    return g;
}


// Exported (module level) Python functions ----------------------------------

static char make_nygraph_docs[] =
"make_nygraph(g): \n\
    Convert the Python NyGraph object into a Nauty/C NyGraph object.\n\
    Return a handle to the Nauty/C NyGraph object as a PyCObject.\n";

static PyObject*
make_nygraph(PyObject *self, PyObject *args)
// Convert the Python NyGraph object into a Nauty/C NyGraph object
// and set Nauty options.
{
    PyObject * py_graph;
    NyGraph * g;

    if (!PyArg_ParseTuple(args, "O", &py_graph)) {
        PyErr_SetString(PyExc_TypeError, "Missing argument.");
        return NULL;
    }

    g = _make_nygraph(py_graph);
    if (g == NULL) return NULL;

    //return PyCObject_FromVoidPtr((void *) g, NULL);
    return PyCapsule_New((void *) g, NULL, NULL);
}


static char delete_nygraph_docs[] =
"delete_nygraph(g): \n\
    Free the allocated memeory for Nauty NyGraph.\n";

static PyObject*
delete_nygraph(PyObject *self, PyObject *args) {
    PyObject *p;
    NyGraph *g;

    if (!PyArg_ParseTuple(args, "O", &p)) {
        PyErr_SetString(PyExc_TypeError, "Missing argument.");
        return NULL;
    }

    //g = (NyGraph *) PyCObject_AsVoidPtr(p);
    g = (NyGraph *) PyCapsule_GetPointer(p, NULL);
    destroy_nygraph(g);
    return Py_BuildValue("");
}


static char graph_autgrp_docs[] =
"graph_autgrp(g):\n\
    Return the (generators, order, orbits, orbit_no)\n\
    of the automorphism group of NyGraph 'g'.\n";

static PyObject*
graph_autgrp(PyObject *self, PyObject *args)
{
    PyObject *py_graph;
    NyGraph * g;
    PyObject *pyret;

    if (!PyArg_ParseTuple(args, "O", &py_graph)) {
        PyErr_SetString(PyExc_TypeError, "Missing argument.");
        return NULL;
    }
    g = _make_nygraph(py_graph);
    if (g == NULL) return NULL;

    // compute automorphism group only
    g->options->getcanon = FALSE;
    g->options->userautomproc = store_generator;

    // initialize storage for collecting generators
    if ((g->generator = malloc(NUM_GENS_INCR * sizeof(permutation*)))
            == NULL) {
        PyErr_SetString(PyExc_MemoryError,
                "Initial generator list allocation failed.");
        return NULL;
    }
    g->max_no_generators = NUM_GENS_INCR;
    
    // *** nauty ***
    // compute automorphism group
    nauty(g->matrix, g->lab, g->ptn, NULL, g->orbits,
            g->options, g->stats,  g->workspace, g->worksize,
            g->no_setwords, g->no_vertices, NULL);
    
    pyret = py_auto_group(g);
    destroy_nygraph(g);
    return pyret;
}


static char graph_cert_docs[] =
"graph_cert(g): \n\
    Return the unique certificate of NyGraph 'g'.\n";

static PyObject*
graph_cert(PyObject *self, PyObject *args)
{
    PyObject *py_graph;
    NyGraph * g;
    PyObject *pyret;

    if (!PyArg_ParseTuple(args, "O", &py_graph)) {
        PyErr_SetString(PyExc_TypeError, "Missing argument.");
        return NULL;
    }
    g = _make_nygraph(py_graph);
    if (g == NULL) return NULL;

    // produce graph certificate by computing canonical labeling
    g->options->getcanon = TRUE;
    if (extend_canonical(g) == NULL) {
        PyErr_SetString(PyExc_MemoryError,
                "Allocating canonical matrix failed");
        return NULL;
    }
    // the produced generators are ignored
    g->options->userautomproc = NULL;

    // *** nauty ***
    nauty(g->matrix, g->lab, g->ptn, NULL, g->orbits,
            g->options, g->stats,  g->workspace, g->worksize,
            g->no_setwords, g->no_vertices, g->cmatrix);

#if PY_MAJOR_VERSION >= 3
    pyret = Py_BuildValue("y#", g->cmatrix,
            g->no_vertices * g->no_setwords * sizeof(setword));
#else
    pyret = Py_BuildValue("s#", g->cmatrix,
            g->no_vertices * g->no_setwords * sizeof(setword));
#endif
    destroy_nygraph(g);
    return pyret;
}


//  Python module initialization  =============================================

static PyMethodDef nautywrap_methods[] = {
    {"graph_cert", graph_cert, METH_VARARGS, graph_cert_docs},
    {"graph_autgrp", graph_autgrp, METH_VARARGS, graph_autgrp_docs},
    {"make_nygraph", make_nygraph, METH_VARARGS, make_nygraph_docs},
    {"delete_nygraph", delete_nygraph, METH_VARARGS, delete_nygraph_docs},
    {NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "nautywrap",
    .m_doc = PyDoc_STR("Graph (auto/iso)morphism wrapper for nauty"),
    .m_size = -1,
    .m_methods = nautywrap_methods,
};

PyObject *
PyInit_nautywrap(void) {
    PyObject *m;

    m = PyModule_Create(&moduledef);
    return m;
#else
void
initnautywrap(void) {
    PyObject *m;

    m = Py_InitModule3("nautywrap", nautywrap_methods,
            "Graph (auto/iso)morphism wrapper for nauty");
#endif
}

