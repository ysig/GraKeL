#include <Python.h>
#include <cstdio>

#include "graph.hh"

static void _destroy(void *g)
{
  if(g)
    {
      //fprintf(stderr, "Free: %x\n", (unsigned int)g);
      delete (bliss::Graph*)g;
    }
}

static PyObject *
graph_create(PyObject *self, PyObject *args)
{
  bliss::Graph *g = new bliss::Graph();
  if(!g)
    Py_RETURN_NONE;
  
  //fprintf(stderr, "Alloc: %x\n", (unsigned int)g);
  PyObject *py_g = PyCObject_FromVoidPtr(g, &_destroy);
  if(!py_g)
    Py_RETURN_NONE;
  return py_g;
}


#if 0
static PyObject *
graph_delete(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;

  if(!PyArg_ParseTuple(args, "O", &py_g))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;
  bliss::Graph *g = (bliss::Graph *)PyCObject_AsVoidPtr(py_g);
  delete g;
  Py_RETURN_NONE;
}
#endif


#if 0
static PyObject *
graph_read_dimacs(PyObject *self, PyObject *args)
{
  const char *filename = 0;

  if(!PyArg_ParseTuple(args, "s", &filename))
    Py_RETURN_NONE;
  FILE *fp = fopen(filename, "r");
  if(!fp)
    Py_RETURN_NONE;
  bliss::Graph *g = bliss::Graph::read_dimacs(fp, 0);
  fclose(fp);
  if(!g)
    Py_RETURN_NONE;
  return PyCObject_FromVoidPtr(g, &_destroy);
}
#endif


#if 0
static PyObject *
pybliss_write_dot(PyObject *self, PyObject *args)
{
  const char *filename = 0;
  PyObject *py_g = NULL;

  if(!PyArg_ParseTuple(args, "Os", &py_g, &filename))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;
  bliss::AbstractGraph *g = (bliss::AbstractGraph *)PyCObject_AsVoidPtr(py_g);
  
  g->write_dot(filename);
  Py_RETURN_NONE;
}
#endif


#if 0
static PyObject *
nof_vertices(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;

  if(!PyArg_ParseTuple(args, "O", &py_g))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;
  bliss::AbstractGraph *g = (bliss::AbstractGraph *)PyCObject_AsVoidPtr(py_g);
  PyObject *py_N = PyInt_FromLong(g->get_nof_vertices());
  if(!py_N)
    Py_RETURN_NONE;
  return py_N;
}
#endif


static PyObject *
add_vertex(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;
  unsigned int color;

  if(!PyArg_ParseTuple(args, "OI", &py_g, &color))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;

  bliss::Graph* g = (bliss::Graph *)PyCObject_AsVoidPtr(py_g);
  assert(g);

  unsigned int v = g->add_vertex(color);
  PyObject *py_N = PyInt_FromLong(v);
  if(!py_N)
    Py_RETURN_NONE;
  return py_N;
}


static PyObject *
add_edge(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;
  unsigned int v1;
  unsigned int v2;

  if(!PyArg_ParseTuple(args, "OII", &py_g, &v1, &v2))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;

  bliss::Graph* g = (bliss::Graph *)PyCObject_AsVoidPtr(py_g);
  assert(g);

  g->add_edge(v1, v2);
  Py_RETURN_NONE;
}


#if 0
static PyObject *
graph_permute(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;
  PyObject *py_perm = NULL;

  /* Read args: graph and permutation objects */
  if(!PyArg_ParseTuple(args, "OO", &py_g, &py_perm))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;
  if(!PyList_Check(py_perm))
    Py_RETURN_NONE;

  bliss::Graph *g = (bliss::Graph *)PyCObject_AsVoidPtr(py_g);
  if(!g)
    Py_RETURN_NONE;

  /* Convert permutation from Python to C */
  const unsigned int N = g->get_nof_vertices();
  if((unsigned int)PyList_Size(py_perm) != N)
    Py_RETURN_NONE;
  unsigned int *perm = (unsigned int*)malloc(N * sizeof(unsigned int));
  if(!perm)
    Py_RETURN_NONE;
  for(unsigned int i = 0; i < N; i++)
    {
      PyObject* py_image = PyList_GetItem(py_perm, i);
      if(!py_image || !PyInt_Check(py_image))
	{
	  free(perm);
	  Py_RETURN_NONE;
	}
      unsigned int image = (unsigned int)PyInt_AsLong(py_image);
      if(image >= N)
	{
	  free(perm);
	  Py_RETURN_NONE;
	}
      perm[i] = image;
    }

  /* Permute the graph */
  bliss::Graph *g2 = g->permute(perm);
  free(perm);
  if(!g2)
    Py_RETURN_NONE;

  return PyCObject_FromVoidPtr(g2, &_destroy);
}
#endif


typedef struct {
  PyObject *py_reporter;
  PyObject *py_reporter_arg;
} ReporterStruct;


static void
_reporter(void *user_param,
	  const unsigned int N,
	  const unsigned int *aut)
{
  if(!user_param)
    return;
  ReporterStruct *s = (ReporterStruct *)user_param;
  /*PyObject *py_reporter = (PyObject*)user_param;*/
  if(!s->py_reporter)
    return;

  PyObject* py_aut = PyList_New(N);
  if(!py_aut)
    {
      return;
    }
  
  for(unsigned int i = 0; i < N; i++)
    if(PyList_SetItem(py_aut, i, PyInt_FromLong((long)aut[i])) != 0)
      return;

  PyObject* args = PyTuple_Pack(2, py_aut, s->py_reporter_arg);
  PyObject* result = PyObject_Call(s->py_reporter, args, NULL);
  if(result)
    Py_DECREF(result);
  Py_DECREF(args);
  Py_DECREF(py_aut);
}



static PyObject *
pybliss_canonical_form(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;
  PyObject *py_reporter = NULL;
  PyObject *py_reporter_arg = NULL;
  bliss::Graph *g = 0;

  if(!PyArg_ParseTuple(args, "OOO", &py_g, &py_reporter, &py_reporter_arg))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;
  if(!PyFunction_Check(py_reporter))
    {
      assert(py_reporter == Py_None);
      py_reporter = NULL;
    }

  g = (bliss::Graph *)PyCObject_AsVoidPtr(py_g);
  assert(g);

  bliss::Stats stats;
  ReporterStruct s;
  s.py_reporter = py_reporter;
  s.py_reporter_arg = py_reporter_arg;
  const unsigned int *cl = g->canonical_form(stats, &_reporter, &s);

  const unsigned int N = g->get_nof_vertices();

  PyObject* py_cl = PyList_New(N);
  if(!py_cl)
    Py_RETURN_NONE;
  
  for(unsigned int i = 0; i < N; i++)
    if(PyList_SetItem(py_cl, i, PyInt_FromLong((long)cl[i])) != 0)
      Py_RETURN_NONE;

  return py_cl;
}


static PyObject *
pybliss_find_automorphisms(PyObject *self, PyObject *args)
{
  PyObject *py_g = NULL;
  PyObject *py_reporter = NULL;
  PyObject *py_reporter_arg = NULL;
  bliss::Graph *g = 0;

  if(!PyArg_ParseTuple(args, "OOO", &py_g, &py_reporter, &py_reporter_arg))
    Py_RETURN_NONE;
  if(!PyCObject_Check(py_g))
    Py_RETURN_NONE;
  if(!PyFunction_Check(py_reporter))
    {
      assert(py_reporter == Py_None);
      py_reporter = NULL;
    }

  g = (bliss::Graph *)PyCObject_AsVoidPtr(py_g);
  assert(g);

  bliss::Stats stats;
  ReporterStruct s;
  s.py_reporter = py_reporter;
  s.py_reporter_arg = py_reporter_arg;
  g->find_automorphisms(stats, &_reporter, &s);
  Py_RETURN_NONE;
}


static PyMethodDef Methods[] = {
    {"create",  graph_create, METH_VARARGS, ""},
    /*{"delete",  graph_delete, METH_VARARGS, ""},*/
    /*{"read_dimacs", graph_read_dimacs, METH_VARARGS, ""},*/
    /*{"nof_vertices", nof_vertices, METH_VARARGS, ""},*/
    {"add_vertex", add_vertex, METH_VARARGS, ""},
    {"add_edge", add_edge, METH_VARARGS, ""},
    /*{"write_dot",  pybliss_write_dot, METH_VARARGS,
      "Write the graph into a file in the graphviz dot format."},*/
    {"canonical_form",  pybliss_canonical_form, METH_VARARGS,
     "Transform the graph into canonical form."},
    {"find_automorphisms",  pybliss_find_automorphisms, METH_VARARGS,
     "Find a generating set for Aut(G)."},
    /*{"graph_permute", graph_permute, METH_VARARGS, ""},*/
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initintpybliss(void)
{
  (void)Py_InitModule("intpybliss", Methods);
}

