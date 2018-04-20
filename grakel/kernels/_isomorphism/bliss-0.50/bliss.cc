#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "defs.hh"
#include "graph.hh"
#include "timer.hh"
#include "utils.hh"

/*
 * Copyright (c) Tommi Junttila
 * Released under the GNU General Public License version 2.
 */

/**
 * \page executable bliss executable
 */

static char *infilename = 0;

static bool opt_directed = false;
static bool opt_canonize = false;
static const char *opt_output_can_file = 0;
static const char *opt_splitting_heuristics = "fm";


static unsigned int verbose_level = 1;
static FILE *verbstr = stdout;


static void usage(FILE *fp, char *argv0)
{
  const char *program_name;
  
  program_name = rindex(argv0, '/');
  
  if(program_name) program_name++;
  else program_name = argv0;
  
  if(!*program_name) program_name = "bliss";
  fprintf(fp, "bliss version %s (compiled "__DATE__")\n", bliss::version);
  fprintf(fp, "Copyright 2003-2007 Tommi Junttila\n");
  fprintf(fp,
"\n"
"Usage: %s [options] [<graph file>]\n"
"\n"
"  -directed   the input graph is directed\n"
"  -can        compute canonical form\n"
"  -ocan=f     compute canonical form and output it in file f\n"
"  -v=N        set verbose level to N; N >= 0 and 1 by default\n"
"  -sh=x       select splitting heuristics, where x is\n"
"                f    first non-singleton cell\n"
"                fl   first largest non-singleton cell\n"
"                fs   first smallest non-singleton cell\n"
"                fm   first maximally non-trivially connected non-singleton cell [default]\n"
"                flm  first largest maximally non-trivially connected non-singleton cell\n"
"                fsm  first smallest maximally non-trivially connected non-singleton cell\n"
"  -version    print the version number and exit\n"
"  -help       print this help and exit\n"
          ,program_name
	  );
}


static void parse_options(int argc, char ** argv)
{
  unsigned int tmp;
  for(int i = 1; i < argc; i++)
    {
      if(strcmp(argv[i], "-can") == 0)
	opt_canonize = true;
      else if((strncmp(argv[i], "-ocan=", 6) == 0) && (strlen(argv[i]) > 6))
	{
	  opt_canonize = true;
	  opt_output_can_file = argv[i]+6;
	}
      else if(sscanf(argv[i], "-v=%u", &tmp) == 1)
	verbose_level = tmp;
      else if(strcmp(argv[i], "-directed") == 0)
	opt_directed = true;
      else if((strncmp(argv[i], "-sh=", 4) == 0) && (strlen(argv[i]) > 4))
	{
	  opt_splitting_heuristics = argv[i]+4;
	}
      else if(strcmp(argv[i], "-version") == 0)
	{
	  fprintf(stdout, "bliss version %s\n", bliss::version);
	  exit(0);
	}
      else if(strcmp(argv[i], "-help") == 0)
	{
	  usage(stdout, argv[0]);
	  exit(0);
	}
      else if(argv[i][0] == '-')
	{
	  fprintf(stderr, "Unknown command line argument `%s'\n", argv[i]);
	  usage(stderr, argv[0]);
	  exit(1);
	}
      else
	{
	  if(infilename)
	    {
	      fprintf(stderr, "Too many file arguments\n");
	      usage(stderr, argv[0]);
	      exit(1);
	    }
	  else
	    {
	      infilename = argv[i];
	    }
	}
    }
}


/**
 * The hook function that prints the found automorphisms.
 * \a param must be a file descriptor (FILE *).
 */
static void report_aut(void *param, unsigned int n, const unsigned int *aut)
{
  assert(param);
  fprintf((FILE*)param, "Generator: ");
  bliss::print_permutation((FILE*)param, n, aut, 1);
  fprintf((FILE*)param, "\n");
}



int main(int argc, char **argv)
{
  bliss::Timer timer;
  bliss::AbstractGraph *g = 0;

  parse_options(argc, argv);

  /* Parse splitting heuristics */
  bliss::Digraph::SplittingHeuristic shs_directed = bliss::Digraph::shs_fm;
  bliss::Graph::SplittingHeuristic shs_undirected = bliss::Graph::shs_fm;
  if(opt_directed)
    {
      if(strcmp(opt_splitting_heuristics, "f") == 0)
	shs_directed = bliss::Digraph::shs_f;
      else if(strcmp(opt_splitting_heuristics, "fs") == 0)
	shs_directed = bliss::Digraph::shs_fs;
      else if(strcmp(opt_splitting_heuristics, "fl") == 0)
	shs_directed = bliss::Digraph::shs_fl;
      else if(strcmp(opt_splitting_heuristics, "fm") == 0)
	shs_directed = bliss::Digraph::shs_fm;
      else if(strcmp(opt_splitting_heuristics, "fsm") == 0)
	shs_directed = bliss::Digraph::shs_fsm;
      else if(strcmp(opt_splitting_heuristics, "flm") == 0)
	shs_directed = bliss::Digraph::shs_flm;
      else
	{
	  fprintf(stderr, "Illegal option -sh=%s, aborting\n",
		  opt_splitting_heuristics);
	  exit(1);
	}
    }
  else
    {
      if(strcmp(opt_splitting_heuristics, "f") == 0)
	shs_undirected = bliss::Graph::shs_f;
      else if(strcmp(opt_splitting_heuristics, "fs") == 0)
	shs_undirected = bliss::Graph::shs_fs;
      else if(strcmp(opt_splitting_heuristics, "fl") == 0)
	shs_undirected = bliss::Graph::shs_fl;
      else if(strcmp(opt_splitting_heuristics, "fm") == 0)
	shs_undirected = bliss::Graph::shs_fm;
      else if(strcmp(opt_splitting_heuristics, "fsm") == 0)
	shs_undirected = bliss::Graph::shs_fsm;
      else if(strcmp(opt_splitting_heuristics, "flm") == 0)
	shs_undirected = bliss::Graph::shs_flm;
      else
	{
	  fprintf(stderr, "Illegal option -sh=%s, aborting\n",
		  opt_splitting_heuristics);
	  exit(1);
	}
    }

  /* Open the input file */
  FILE *infile = stdin;
  if(infilename)
    {
      infile = fopen(infilename, "r");
      if(!infile)
	{
	  fprintf(stderr, "Cannot not open `%s' for input, aborting\n",
		  infilename);
	  exit(1);
	}
    }

  /* Read the graph from the file */
  if(opt_directed)
    {
      /* Read directed graph in the DIMACS format */
      g = bliss::Digraph::read_dimacs(infile);
    }
  else
    {
      /* Read undirected graph in the DIMACS format */
      g = bliss::Graph::read_dimacs(infile);
    }
  
  if(infile != stdin)
    fclose(infile);

  if(!g)
    {
      /* Failed to read the graph, exit with non-zero status. */
      exit(1);
    }
  
  if(verbose_level >= 2)
    {
      fprintf(verbstr, "Graph read in %.2f seconds\n", timer.get_duration());
      fflush(verbstr);
    }


  bliss::Stats stats;

  /* Set splitting heuristics and verbose level */
  if(opt_directed)
    ((bliss::Digraph*)g)->set_splitting_heuristic(shs_directed);
  else
    ((bliss::Graph*)g)->set_splitting_heuristic(shs_undirected);
  g->set_verbose_level(verbose_level);
  g->set_verbose_file(verbstr);


  if(opt_canonize)
    {
      const unsigned int *cl = g->canonical_form(stats, &report_aut, stdout);

      fprintf(stdout, "Canonical labeling: ");
      bliss::print_permutation(stdout, g->get_nof_vertices(), cl, 1);
      fprintf(stdout, "\n");

      if(opt_output_can_file)
	{
	  bliss::AbstractGraph *cf = g->permute(cl);
	  FILE *fp = fopen(opt_output_can_file, "w");
	  if(!fp)
	    {
	      fprintf(stderr,
		      "Cannot open '%s' for outputting the canonical form",
		      opt_output_can_file);
	      exit(1);
	    }
	  cf->write_dimacs(fp);
	  fclose(fp);
	  delete cf;
	}
    }
  else
    {
      g->find_automorphisms(stats, &report_aut, stdout);
    }

  if(verbose_level > 0)
    {
      fprintf(verbstr, "Nodes:\t\t%lu\n", stats.nof_nodes);
      fprintf(verbstr, "Leaf nodes:\t%lu\n", stats.nof_leaf_nodes);
      fprintf(verbstr, "Bad nodes:\t%lu\n", stats.nof_bad_nodes);
      fprintf(verbstr, "Canrep updates:\t%lu\n", stats.nof_canupdates);
      fprintf(verbstr, "Generators:\t%lu\n", stats.nof_generators);
      fprintf(verbstr, "Max level:\t%lu\n", stats.max_level);
      fprintf(verbstr, "|Aut|:\t\t");
      stats.group_size.print(verbstr);
      fprintf(verbstr, "\n");
      fflush(verbstr);
    }

  if(verbose_level > 0)
    {
      fprintf(verbstr, "Total time:\t%.2f seconds\n", timer.get_duration());
      fflush(verbstr);
    }

  delete g; g = 0;

  return 0;
}
