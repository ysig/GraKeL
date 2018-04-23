#include <cstdio>
#include <cassert>
#include <climits>
#include <set>
#include <list>
#include <algorithm>
#include "defs.hh"
#include "graph.hh"
#include "partition.hh"
#include "utils.hh"

/*
 * Copyright (c) Tommi Junttila
 * Released under the GNU General Public License version 2.
 */

namespace bliss {

static const bool not_yet_implemented = false;
static const bool should_not_happen = false;

/*-------------------------------------------------------------------------
 *
 * Constructor and destructor routines for the abstract graph class
 *
 *-------------------------------------------------------------------------*/


AbstractGraph::AbstractGraph()
{
  /* Initialize stuff */
  first_path_labeling = 0;
  first_path_labeling_inv = 0;
  best_path_labeling = 0;
  best_path_labeling_inv = 0;
  first_path_automorphism = 0;
  best_path_automorphism = 0;
  in_search = false;

  opt_use_long_prune = true;


  verbose_level = 0;
  verbstr = stdout;

  report_hook = 0;
  report_user_param = 0;
}


AbstractGraph::~AbstractGraph()
{
  if(first_path_labeling) {
    free(first_path_labeling); first_path_labeling = 0; }
  if(first_path_labeling_inv) {
    free(first_path_labeling_inv); first_path_labeling_inv = 0; }
  if(best_path_labeling) {
    free(best_path_labeling); best_path_labeling = 0; }
  if(best_path_labeling_inv) {
    free(best_path_labeling_inv); best_path_labeling_inv = 0; }
  if(first_path_automorphism) {
    free(first_path_automorphism); first_path_automorphism = 0; }
  if(best_path_automorphism) {
    free(best_path_automorphism); best_path_automorphism = 0; }
 
  report_hook = 0;
  report_user_param = 0;
}

void AbstractGraph::set_verbose_level(const unsigned int level)
{
  verbose_level = level;
}

void AbstractGraph::set_verbose_file(FILE * const fp)
{
  verbstr = fp;
}



/*-------------------------------------------------------------------------
 *
 * Routines for refinement to equitable partition
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::refine_to_equitable()
{
  BLISS_ASSERT(p.splitting_queue.is_empty());

  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      p.add_in_splitting_queue(cell);
    }

  return do_refine_to_equitable();
}


void AbstractGraph::refine_to_equitable(Partition::Cell * const unit_cell)
{
  BLISS_ASSERT(unit_cell->is_unit());

#if defined(BLISS_EXPENSIVE_CONSISTENCY_CHECKS)
  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      assert(cell->in_splitting_queue == false);
      assert(cell->in_neighbour_heap == false);
    }
#endif
  
  BLISS_ASSERT(p.splitting_queue.is_empty());

  p.add_in_splitting_queue(unit_cell);

  do_refine_to_equitable();
  p.consistency_check();

}



void AbstractGraph::refine_to_equitable(Partition::Cell * const unit_cell1,
					Partition::Cell * const unit_cell2)
{
  BLISS_ASSERT(unit_cell1->is_unit());
  BLISS_ASSERT(unit_cell2->is_unit());

#if defined(BLISS_EXPENSIVE_CONSISTENCY_CHECKS)
  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      assert(cell->in_splitting_queue == false);
      assert(cell->in_neighbour_heap == false);
    }
#endif
  
  BLISS_ASSERT(p.splitting_queue.is_empty());

  p.add_in_splitting_queue(unit_cell1);
  p.add_in_splitting_queue(unit_cell2);

  do_refine_to_equitable();
  p.consistency_check();

}



void AbstractGraph::do_refine_to_equitable()
{
  assert(!p.splitting_queue.is_empty());
  assert(neighbour_heap.is_empty());

  eqref_hash.reset();

  while(!p.splitting_queue.is_empty())
    {
      Partition::Cell *cell = p.splitting_queue.pop_front();
      BLISS_ASSERT(cell->in_splitting_queue);
      cell->in_splitting_queue = false;

      if(cell->is_unit())
	{
	  if(in_search) {
	    if(first_path_automorphism) {
	      /* Build the (potential) automorphism on-the-fly */
	      BLISS_ASSERT(first_path_labeling_inv);
	      first_path_automorphism[first_path_labeling_inv[cell->first]] =
		p.elements[cell->first];
	    }
	    if(best_path_automorphism)
	      {
		/* Build the (potential) automorphism on-the-fly */
		BLISS_ASSERT(best_path_labeling_inv);
		best_path_automorphism[best_path_labeling_inv[cell->first]] =
		  p.elements[cell->first];
	      }
	  }
	  
	  bool worse = split_neighbourhood_of_unit_cell(cell);
	  if(in_search && worse)
	    goto worse_exit;
	}
      else
	{
	  split_neighbourhood_of_cell(cell);
	}
    }

  eqref_worse_than_certificate = false;
  return;

 worse_exit:
  /* Clear splitting_queue */
  p.clear_splitting_queue();
  eqref_worse_than_certificate = true;
  return;
}





/*-------------------------------------------------------------------------
 *
 * Routines for handling the canonical labeling
 *
 *-------------------------------------------------------------------------*/

/** \internal
 * Assign the labeling induced by the current partition 'this.p' to
 * \a labeling.
 * That is, if the partition is [[2,0],[1]],
 * then \a labeling will map 0 to 1, 1 to 2, and 2 to 0.
 */
void AbstractGraph::update_labeling(unsigned int * const labeling)
{
  const unsigned int N = get_nof_vertices();
  unsigned int *ep = p.elements;
  for(unsigned int i = 0; i < N; i++, ep++)
    labeling[*ep] = i;
}

/** \internal
 * The same as update_labeling() except that the inverse of the labeling
 * is also produced and assigned to \a labeling_inv.
 */
void AbstractGraph::update_labeling_and_its_inverse(unsigned int * const labeling,
						    unsigned int * const labeling_inv)
{
  const unsigned int N = get_nof_vertices();
  unsigned int *ep = p.elements;
  unsigned int *clip = labeling_inv;

  for(unsigned int i = 0; i < N; ) {
    labeling[*ep] = i;
    i++;
    *clip = *ep;
    ep++;
    clip++;
  }
}





/*-------------------------------------------------------------------------
 *
 * Routines for handling automorphisms
 *
 *-------------------------------------------------------------------------*/


/** \internal
 * Reset the permutation \a perm to the identity permutation.
 */
void AbstractGraph::reset_permutation(unsigned int *perm)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++, perm++)
    *perm = i;
}


bool AbstractGraph::is_automorphism(unsigned int * const perm)
{
  assert(should_not_happen);
  return false;
}





/*-------------------------------------------------------------------------
 *
 * Long prune code
 *
 *-------------------------------------------------------------------------*/

void AbstractGraph::long_prune_init()
{
  const unsigned int N = get_nof_vertices();
  long_prune_temp.clear();
  long_prune_temp.resize(N);
#if defined(BLISS_CONSISTENCY_CHECKS)
  for(unsigned int i = 0; i < N; i++)
    assert(long_prune_temp[i] == false);
#endif
  /* Of how many automorphisms we can store information in
     the predefined, fixed amount of memory? */
  const unsigned int nof_fitting_in_max_mem =
    (long_prune_options_max_mem * 1024 * 1024) / (((N * 2) / 8)+1);
  long_prune_max_stored_autss = long_prune_options_max_stored_auts;
  /* Had some problems with g++ in using (a<b)?a:b when constants involved,
     so had to make this in a stupid way... */
  if(nof_fitting_in_max_mem < long_prune_options_max_stored_auts)
    long_prune_max_stored_autss = nof_fitting_in_max_mem;

  while(!long_prune_fixed.empty())
    {
      delete long_prune_fixed.back();
      long_prune_fixed.pop_back();
    }
  while(!long_prune_mcrs.empty())
    {
      delete long_prune_mcrs.back();
      long_prune_mcrs.pop_back();
    }
  long_prune_fixed.resize(N, 0);
  long_prune_mcrs.resize(N, 0);
#if 0
  for(unsigned int i = 0; i < long_prune_max_stored_autss; i++)
    {
      long_prune_fixed.push_back(new std::vector<bool>(N));
      long_prune_mcrs.push_back(new std::vector<bool>(N));
    }
#endif
  long_prune_begin = 0;
  long_prune_end = 0;
}


void AbstractGraph::long_prune_deallocate()
{
  while(!long_prune_fixed.empty())
    {
      delete long_prune_fixed.back();
      long_prune_fixed.pop_back();
    }
  while(!long_prune_mcrs.empty())
    {
      delete long_prune_mcrs.back();
      long_prune_mcrs.pop_back();
    }
}


void AbstractGraph::long_prune_swap(const unsigned int i, const unsigned int j)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(long_prune_begin <= i);
  assert(i < long_prune_end);
  assert(long_prune_begin <= j);
  assert(j < long_prune_end);
  const unsigned int real_i = i % long_prune_max_stored_autss;
  const unsigned int real_j = j % long_prune_max_stored_autss;
  assert(long_prune_fixed[real_i]);
  assert(long_prune_fixed[real_j]);
  assert(long_prune_mcrs[real_i]);
  assert(long_prune_mcrs[real_j]);
  std::vector<bool> * tmp = long_prune_fixed[real_i];
  long_prune_fixed[real_i] = long_prune_fixed[real_j];
  long_prune_fixed[real_j] = tmp;
  tmp = long_prune_mcrs[real_i];
  long_prune_mcrs[real_i] = long_prune_mcrs[real_j];
  long_prune_mcrs[real_j] = tmp;
}


std::vector<bool> &
AbstractGraph::long_prune_allocget_fixed(const unsigned int index)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(long_prune_begin <= index);
  assert(index < long_prune_end);
  const unsigned int i = index % long_prune_max_stored_autss;
  if(!long_prune_fixed[i])
    long_prune_fixed[i] = new std::vector<bool>(get_nof_vertices());
  return *long_prune_fixed[i];
}


std::vector<bool> &
AbstractGraph::long_prune_get_fixed(const unsigned int index)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(long_prune_begin <= index);
  assert(index < long_prune_end);
  return *long_prune_fixed[index % long_prune_max_stored_autss];
}


std::vector<bool> &
AbstractGraph::long_prune_allocget_mcrs(const unsigned int index)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(long_prune_begin <= index);
  assert(index < long_prune_end);
  const unsigned int i = index % long_prune_max_stored_autss;
  if(!long_prune_mcrs[i])
    long_prune_mcrs[i] = new std::vector<bool>(get_nof_vertices());
  return *long_prune_mcrs[i];
}


std::vector<bool> &
AbstractGraph::long_prune_get_mcrs(const unsigned int index)
{
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  assert(long_prune_begin <= index);
  assert(index < long_prune_end);
  return *long_prune_mcrs[index % long_prune_max_stored_autss];
}


void AbstractGraph::long_prune_add_automorphism(const unsigned int *aut)
{
  if(long_prune_max_stored_autss == 0)
    return;

  const unsigned int N = get_nof_vertices();

#if defined(BLISS_CONSISTENCY_CHECKS)
  assert(long_prune_temp.size() == N);
  for(unsigned int i = 0; i < N; i++)
    assert(long_prune_temp[i] == false);
#endif

  BLISS_ASSERT(long_prune_fixed.size() == long_prune_mcrs.size());
  assert(long_prune_begin <= long_prune_end);
  assert(long_prune_end - long_prune_end <= long_prune_max_stored_autss);
  if(long_prune_end - long_prune_begin == long_prune_max_stored_autss)
    {
      long_prune_begin++;
    }
  long_prune_end++;
  std::vector<bool> &fixed = long_prune_allocget_fixed(long_prune_end-1);
  std::vector<bool> &mcrs = long_prune_allocget_mcrs(long_prune_end-1);
  /* Mark nodes that are (i) fixed or (ii) minimal orbit representatives
   * under the automorphism 'aut' */
  for(unsigned int i = 0; i < N; i++)
    {
      fixed[i] = (aut[i] == i);
      if(!long_prune_temp[i])
	{
	  mcrs[i] = true;
	  unsigned int j = aut[i];
	  while(j != i)
	    {
	      assert(i <= j);
	      long_prune_temp[j] = true;
	      j = aut[j];
	    }
	}
      else
	{
	  mcrs[i] = false;
	}
      long_prune_temp[i] = false;
    }


#if defined(BLISS_CONSISTENCY_CHECKS)
  for(unsigned int i = 0; i < N; i++)
    assert(long_prune_temp[i] == false);
#endif
}


/*-------------------------------------------------------------------------
 *
 * Routines for handling orbit information
 *
 *-------------------------------------------------------------------------*/


void AbstractGraph::update_orbit_information(Orbit &o, const unsigned int *p)
{
  const unsigned int N = get_nof_vertices();
  for(unsigned int i = 0; i < N; i++)
    if(p[i] != i)
      o.merge_orbits(i, p[i]);
}






/*-------------------------------------------------------------------------
 *
 * The actual backtracking search
 *
 *-------------------------------------------------------------------------*/


typedef struct
{
  int split_element;
  unsigned int split_cell_first;
  Partition::BacktrackPoint partition_bt_point;
  unsigned int certificate_index;

  bool in_first_path;
  bool in_best_path;
  bool equal_to_first_path;
  int cmp_to_best_path;

  bool needs_long_prune;
  unsigned int long_prune_begin;
  std::set<unsigned int, std::less<unsigned int> > long_prune_redundant;
  
  UintSeqHash eqref_hash;
  unsigned int subcertificate_length;
} LevelInfo;



typedef struct {
  unsigned int splitting_element;
  unsigned int certificate_index;
  unsigned int subcertificate_length;
  UintSeqHash eqref_hash;
} PathInfo;


void AbstractGraph::search(const bool canonical,
			   Stats &stats)
{
  const unsigned int N = get_nof_vertices();

  unsigned int all_same_level = UINT_MAX;

  p.graph = this;

  /*
   * Must be done!
   */
  remove_duplicate_edges();

  /*
   * Reset search statistics
   */
  stats.group_size.assign(1);
  stats.group_size_approx = 1.0;
  stats.nof_nodes = 1;
  stats.nof_leaf_nodes = 1;
  stats.nof_bad_nodes = 0;
  stats.nof_canupdates = 0;
  stats.nof_generators = 0;
  stats.max_level = 0;

  if(first_path_labeling)
    {
      free(first_path_labeling);
      first_path_labeling = 0;
    }
  if(first_path_labeling_inv)
    {
      free(first_path_labeling_inv);
      first_path_labeling_inv = 0;
   }
  if(first_path_automorphism)
    {
      free(first_path_automorphism);
      first_path_automorphism = 0;
   }

  if(best_path_labeling)
    {
      free(best_path_labeling);
      best_path_labeling = 0;
    }
  if(best_path_labeling_inv)
    {
      free(best_path_labeling_inv);
      best_path_labeling_inv = 0;
   }
  if(best_path_automorphism)
    {
      free(best_path_automorphism);
      best_path_automorphism = 0;
   }

  if(N == 0)
    return;

  p.init(N);
  neighbour_heap.init(N);

  in_search = false;
  /* eqref_hash is not exploited when building the initial partition,
     do not compute it as it seems to take some time */
  compute_eqref_hash = false;

  p.level = 0;

  make_initial_equitable_partition();

#if defined(BLISS_VERIFY_EQUITABLEDNESS)
  assert(is_equitable() && "Internal error: initial partition not equitable");
#endif

  if(verbstr && verbose_level >= 2)
    {
      //fprintf(verbstr, "Initial partition computed in %.2f seconds\n"); 
      // TIMER.HH removed for allowing windows
      fflush(verbstr);
      //p.print(verbstr); fprintf(verbstr, "\n");
      //write_dot("DEBUG.dot");
    }
  


  /*
   * Allocate space for the labelings
   */
  if(first_path_labeling)
    free(first_path_labeling);
  first_path_labeling = (unsigned int*)calloc(N, sizeof(unsigned int));
  if(best_path_labeling)
    free(best_path_labeling);
  best_path_labeling = (unsigned int*)calloc(N, sizeof(unsigned int));

  /*
   * Are there any non-singleton cells?
   */
  if(p.is_discrete())
    {
      update_labeling(best_path_labeling);
      return;
    }


  /*
   * Allocate space for the inverses of the labelings
   */
  if(first_path_labeling_inv)
    free(first_path_labeling_inv);
  first_path_labeling_inv = (unsigned int*)calloc(N, sizeof(unsigned int));
  if(best_path_labeling_inv)
    free(best_path_labeling_inv);
  best_path_labeling_inv = (unsigned int*)calloc(N, sizeof(unsigned int));


  /*
   * Allocate space for the automorphisms
   */
  if(first_path_automorphism) free(first_path_automorphism);
  first_path_automorphism = (unsigned int*)malloc(N * sizeof(unsigned int));
  if(best_path_automorphism) free(best_path_automorphism);
  best_path_automorphism = (unsigned int*)malloc(N * sizeof(unsigned int));


  /*
   * Initialize orbit information
   */
  first_path_orbits.init(N);
  best_path_orbits.init(N);

  /*
   * Initialize certificate memory
   */
  initialize_certificate();
  //assert(certificate);
  assert(certificate_index == 0);

  LevelInfo info;
  std::vector<LevelInfo> search_stack;
  std::vector<PathInfo> first_path_info;
  std::vector<PathInfo> best_path_info;

  search_stack.clear();
  assert(neighbour_heap.is_empty());

  if(opt_use_long_prune)
    {
      /* Initialize long prune */
      long_prune_init();
    }


  /*
   * Build the first level info
   */
  info.split_cell_first = find_next_cell_to_be_splitted(p.first_cell)->first;
  info.split_element = -1;
  info.partition_bt_point = p.set_backtrack_point();
  info.certificate_index = 0;
  info.in_first_path = false;
  info.in_best_path = false;
  info.long_prune_begin = 0;
  search_stack.push_back(info);

  /*
   * Set status and global flags for search related procedures
   */
  in_search = true;
  compute_eqref_hash = true;
  refine_compare_certificate = false;
  stats.nof_leaf_nodes = 0;



  p.consistency_check();

  /*
   * The actual backtracking search
   */
  while(!search_stack.empty()) 
    {
      info = search_stack.back();
      search_stack.pop_back();

      p.consistency_check();

      /* p.print(stderr); fprintf(stderr, "\n"); */
      /*
       * Restore partition, certificate index, and split cell
       */
      p.goto_backtrack_point(info.partition_bt_point);
      info.partition_bt_point = p.set_backtrack_point();
      assert(info.certificate_index <= certificate_size);
      certificate_index = info.certificate_index;
      refine_current_path_certificate_index = info.certificate_index;
      //certificate_current_path.resize(certificate_index);
      Partition::Cell * const cell = p.element_to_cell_map[p.elements[info.split_cell_first]];
      assert(cell->length > 1);

      p.consistency_check();



      if(p.level > 0 && !info.in_first_path)
	{
	  unsigned int begin = (info.long_prune_begin>long_prune_begin)?info.long_prune_begin:long_prune_begin;
	  for(unsigned int i = begin; i < long_prune_end; i++)
	    {
	      const std::vector<bool> &fixed = long_prune_get_fixed(i);
#if defined(BLISS_CONSISTENCY_CHECKS)
	      for(unsigned int l = 0; l < p.level - 1; l++)
		assert(fixed[search_stack[l].split_element]);
#endif
	      if(fixed[search_stack[p.level-1].split_element] == false)
		{
		  long_prune_swap(begin, i);
		  begin++;
		  info.long_prune_begin = begin;
		  continue;
		}
	    }

	  if(info.split_element == -1)
	    {
	      info.needs_long_prune = true;
	    }
	  else if(info.needs_long_prune)
	    {
	      info.needs_long_prune = false;
	      /* THIS IS A QUITE HORRIBLE HACK! */
	      unsigned int begin = (info.long_prune_begin>long_prune_begin)?info.long_prune_begin:long_prune_begin;
	      for(unsigned int i = begin; i < long_prune_end; i++)
		{
		  const std::vector<bool> &fixed = long_prune_get_fixed(i);
#if defined(BLISS_CONSISTENCY_CHECKS)
		  for(unsigned int l = 0; l < p.level-1; l++)
		    assert(fixed[search_stack[l].split_element]);
#endif
		  if(fixed[search_stack[p.level-1].split_element] == false)
		    {
		      long_prune_swap(begin, i);
		      begin++;
		      info.long_prune_begin = begin;
		      continue;
		    }
		  const std::vector<bool> &mcrs = long_prune_get_mcrs(i);
		  unsigned int *ep = p.elements + cell->first;
		  for(unsigned int j = cell->length; j > 0; j--, ep++) {
		    if(mcrs[*ep] == false)
		      {
			info.long_prune_redundant.insert(*ep);
		      }
		  }
		}
	    }
	}

      /*
       * Find the next smallest element in cell
       */
      unsigned int next_split_element = UINT_MAX;
      unsigned int *next_split_element_pos = 0;
      unsigned int *ep = p.elements + cell->first;
      if(info.in_first_path)
	{
	  /* Find the next larger splitting element that is a mor */
	  for(unsigned int i = cell->length; i > 0; i--, ep++) {
	    if((int)(*ep) > info.split_element &&
	       *ep < next_split_element &&
	       first_path_orbits.is_minimal_representative(*ep)) {
	      next_split_element = *ep;
	      next_split_element_pos = ep;
	    }
	  }
	}
      else if(info.in_best_path)
	{
	  /* Find the next larger splitting element that is a mor */
	  for(unsigned int i = cell->length; i > 0; i--, ep++) {
	    if((int)(*ep) > info.split_element &&
	       *ep < next_split_element &&
	       best_path_orbits.is_minimal_representative(*ep) &&
	       (!opt_use_long_prune ||
		info.long_prune_redundant.find(*ep) ==
		info.long_prune_redundant.end())) {
	      next_split_element = *ep;
	      next_split_element_pos = ep;
	    }
	  }
	}
      else
	{
	  /* Find the next larger splitting element */
	  for(unsigned int i = cell->length; i > 0; i--, ep++) {
	    if((int)(*ep) > info.split_element &&
	       *ep < next_split_element &&
	       (!opt_use_long_prune ||
		info.long_prune_redundant.find(*ep) ==
		info.long_prune_redundant.end())) {
	      next_split_element = *ep;
	      next_split_element_pos = ep;
	    }
	  }
	}
      if(next_split_element == UINT_MAX)
	{
	  /*
	   * No more splitting elements (unexplored children) in the cell
	   */
	  /* Update group size if required */
	  if(info.in_first_path == true) {
	    const unsigned int index =
	      first_path_orbits.orbit_size(first_path_info[p.level].splitting_element);
	    stats.group_size.multiply(index);
	    stats.group_size_approx *= (long double)index;
	    /*
	     * Update all_same_level
	     */
	    if(index == cell->length && all_same_level == p.level+1)
	      all_same_level = p.level;
	    if(verbstr && verbose_level >= 2)
	      {
		fprintf(verbstr,
			"Level %u: orbits=%u, index=%u/%u, all_same_level=%u\n",
			p.level,
			first_path_orbits.nof_orbits(),
			index, cell->length,
			all_same_level);
		fflush(stdout);
	      }
	  }
	  /* Backtrack to the previous level */
	  p.level--;
	  continue;
	}

      /* Split on smallest */
      info.split_element = next_split_element;
      
      /*
       * Save the current search situation
       */
      search_stack.push_back(info);

      /*
       * No more in the first path
       */
      info.in_first_path = false;
      /*
       * No more in the best path
       */
      info.in_best_path = false;

      p.level++;
      stats.nof_nodes++;
      if(p.level > stats.max_level)
	stats.max_level = p.level;

      p.consistency_check();

      /* Individualize, i.e. split the cell in two, the latter new cell
       * will be a unit one containing next_split_element */
      Partition::Cell * const new_cell = p.individualize(cell,
							 next_split_element);


      if(!first_path_info.empty())
	{
	  refine_equal_to_first = info.equal_to_first_path;
	  if(refine_equal_to_first)
	    refine_first_path_subcertificate_end =
	      first_path_info[p.level-1].certificate_index +
	      first_path_info[p.level-1].subcertificate_length;
	  if(canonical)
	    {
	      refine_cmp_to_best = info.cmp_to_best_path;
	      if(refine_cmp_to_best == 0)
		refine_best_path_subcertificate_end =
		  best_path_info[p.level-1].certificate_index +
		  best_path_info[p.level-1].subcertificate_length;
	    }
	  else
	    refine_cmp_to_best = -1;
	}
      /*
       * Refine the new partition to equitable
       */
      if(cell->is_unit())
	refine_to_equitable(cell, new_cell);
      else 
	refine_to_equitable(new_cell);

      p.consistency_check();



      if(p.is_discrete())
	{
	  /* Update statistics */
	  stats.nof_leaf_nodes++;
	}

      if(!first_path_info.empty())
	{
	  /* We are no longer on the first path */
	  assert(best_path_info.size() > 0);
	  assert(refine_current_path_certificate_index >= certificate_index);
	  const unsigned int subcertificate_length = 
	    refine_current_path_certificate_index - certificate_index;
	  if(refine_equal_to_first)
	    {
	      /* Was equal to the first path so far */
	      assert(first_path_info.size() >= p.level);
	      PathInfo &first_pinfo = first_path_info[p.level-1];
	      assert(first_pinfo.certificate_index == certificate_index);
	      if(subcertificate_length != first_pinfo.subcertificate_length)
		{
		  refine_equal_to_first = false;
		}
	      else if(first_pinfo.eqref_hash.cmp(eqref_hash) != 0)
		{
		  refine_equal_to_first = false;
		}
	    }
	  if(canonical && (refine_cmp_to_best == 0))
	    {
	      /* Was equal to the best path so far */
	      assert(best_path_info.size() >= p.level);
	      PathInfo &best_pinfo = best_path_info[p.level-1];
	      assert(best_pinfo.certificate_index == certificate_index);
	      if(subcertificate_length < best_pinfo.subcertificate_length)
		{
		  refine_cmp_to_best = -1;
		}
	      else if(subcertificate_length > best_pinfo.subcertificate_length)
		{
		  refine_cmp_to_best = 1;
		}
	      else if(best_pinfo.eqref_hash.cmp(eqref_hash) > 0)
		{
		  refine_cmp_to_best = -1;
		}
	      else if(best_pinfo.eqref_hash.cmp(eqref_hash) < 0)
		{
		  refine_cmp_to_best = 1;
		}
	    }


	  if(refine_equal_to_first == false &&
	     (!canonical || (refine_cmp_to_best < 0)))
	    {
	      /* Backtrack */
	      stats.nof_bad_nodes++;
	      if(search_stack.back().equal_to_first_path == true &&
		 p.level > all_same_level)
		{
		  assert(all_same_level >= 1);
		  for(unsigned int i = all_same_level;
		      i < search_stack.size();
		      i++)
		    {
		      search_stack[i].equal_to_first_path = false;
		    }
		}


	      while(!search_stack.empty())
		{
		  p.level--;
		  LevelInfo &info2 = search_stack.back();
		  if(!(info2.equal_to_first_path == false &&
		       (!canonical || (info2.cmp_to_best_path < 0))))
		    break;
		  search_stack.pop_back();
		}
	      continue;
	    }
	}

#if defined(BLISS_VERIFY_EQUITABLEDNESS)
      /* The new partition should be equitable */
      assert(is_equitable() &&
	     "Internal error: partition after refinement was not equitable");
#endif

      info.equal_to_first_path = refine_equal_to_first;
      info.cmp_to_best_path = refine_cmp_to_best;

      certificate_index = refine_current_path_certificate_index;

      search_stack.back().eqref_hash = eqref_hash;
      search_stack.back().subcertificate_length =
	certificate_index - info.certificate_index;


      if(!p.is_discrete())
	{
	  /*
	   * An internal, non-leaf node
	   */
	  /* Build the next node info */
	  /* Find the next cell to be splitted */
	  Partition::Cell * const next_split_cell = find_next_cell_to_be_splitted(p.element_to_cell_map[p.elements[info.split_cell_first]]);
	  assert(next_split_cell);
	  /* Copy current info to the search stack */
	  search_stack.push_back(info);
	  LevelInfo &new_info = search_stack.back();
	  new_info.split_cell_first = next_split_cell->first;
	  new_info.split_element = -1;
	  new_info.certificate_index = certificate_index;
	  new_info.partition_bt_point = p.set_backtrack_point();
	  new_info.long_prune_redundant.clear();
	  new_info.long_prune_begin = info.long_prune_begin;
	  continue;
	}

      /*
       * A leaf node
       */
      assert(certificate_index == certificate_size);

      if(first_path_info.empty())
	{
	  /* The first path, update first_path and best_path */
	  //fprintf(stdout, "Level %u: FIRST\n", p.level); fflush(stdout);
	  stats.nof_canupdates++;
	  /*
	   * Update labelings and their inverses
	   */
	  update_labeling_and_its_inverse(first_path_labeling,
					  first_path_labeling_inv);
	  update_labeling_and_its_inverse(best_path_labeling,
					  best_path_labeling_inv);
	  /*
	   * Reset automorphism array
	   */
	  reset_permutation(first_path_automorphism);
	  reset_permutation(best_path_automorphism);
	  /*
	   * Reset orbit information
	   */
	  first_path_orbits.reset();
	  best_path_orbits.reset();
	  /*
	   * Reset group size
	   */
	  stats.group_size.assign(1);
	  stats.group_size_approx = 1.0;
	  /*
	   * Reset all_same_level
	   */
	  all_same_level = p.level;
	  /*
	   * Mark the current path to be the first and best one and save it
	   */
	  const unsigned int base_size = search_stack.size();
	  assert(p.level == base_size);
	  best_path_info.clear();
	  //fprintf(stdout, " New base is: ");
	  for(unsigned int i = 0; i < base_size; i++) {
	    search_stack[i].in_first_path = true;
	    search_stack[i].in_best_path = true;
	    search_stack[i].equal_to_first_path = true;
	    search_stack[i].cmp_to_best_path = 0;
	    PathInfo path_info;
	    path_info.splitting_element = search_stack[i].split_element;
	    path_info.certificate_index = search_stack[i].certificate_index;
	    path_info.eqref_hash = search_stack[i].eqref_hash;
	    path_info.subcertificate_length = search_stack[i].subcertificate_length;
	    first_path_info.push_back(path_info);
	    best_path_info.push_back(path_info);
	    //fprintf(stdout, "%u ", search_stack[i].split_element);
	  }
	  //fprintf(stdout, "\n"); fflush(stdout);
	  certificate_first_path = certificate_current_path;
	  certificate_best_path = certificate_current_path;

	  refine_compare_certificate = true;


	  /*
	   * Backtrack to the previous level
	   */
	  p.level--;
	  continue;
	}

      BLISS_ASSERT(first_path_info.size() > 0);

      //fprintf(stdout, "Level %u: LEAF %d %d\n", p.level, info.equal_to_first_path, info.cmp_to_best_path); fflush(stdout);

      if(info.equal_to_first_path)
	{
	  /*
	   * An automorphism found: aut[i] = elements[first_path_labeling[i]]
	   */
	  assert(!info.in_first_path);
	  //fprintf(stdout, "A"); fflush(stdout);
	  
	  
#if defined(BLISS_CONSISTENCY_CHECKS)
	  /* Verify that the automorphism is correctly built */
	  for(unsigned int i = 0; i < N; i++)
	    assert(first_path_automorphism[i] ==
		   p.elements[first_path_labeling[i]]);
#endif
	  
#if defined(BLISS_VERIFY_AUTOMORPHISMS)
	  /* Verify that it really is an automorphism */
	  assert(is_automorphism(first_path_automorphism));
#endif

	  if(opt_use_long_prune)
	    {
	      long_prune_add_automorphism(first_path_automorphism);
	    }

	  /*
	   * Update orbit information
	   */
	  update_orbit_information(first_path_orbits, first_path_automorphism);
	  
	  /*
	   * Compute backjumping level
	   */
	  unsigned int backjumping_level = 0;
	  for(unsigned int i = search_stack.size(); i > 0; i--) {
	    const unsigned int split_element =
	      search_stack[backjumping_level].split_element;
	    if(first_path_automorphism[split_element] != split_element)
	      break;
	    backjumping_level++;
	  }
	  assert(backjumping_level < p.level);
	  /*
	   * Go back to backjumping_level
	   */
	  p.level = backjumping_level;
	  search_stack.resize(p.level + 1);
	  
	  /* Report automorphism */
	  if(report_hook)
	    (*report_hook)(report_user_param,
			   get_nof_vertices(),
			   first_path_automorphism);
	  /* Update statistics */
	  stats.nof_generators++;
	  continue;
	}

      assert(canonical);
      assert(info.cmp_to_best_path >= 0);
      if(info.cmp_to_best_path > 0)
	{
	  /*
	   * A new, better representative found
	   */
	  //fprintf(stdout, "Level %u: NEW BEST\n", p.level); fflush(stdout);
	  stats.nof_canupdates++;
	  /*
	   * Update canonical labeling and its inverse
	   */
	  update_labeling_and_its_inverse(best_path_labeling,
					  best_path_labeling_inv);
	  /* Reset best path automorphism */
	  reset_permutation(best_path_automorphism);
	  /* Reset best path orbit structure */
	  best_path_orbits.reset();
	  /*
	   * Mark the current path to be the best one and save it
	   */
	  const unsigned int base_size = search_stack.size();
	  assert(p.level == base_size);
	  best_path_info.clear();
	  //fprintf(stdout, " New base is: ");
	  for(unsigned int i = 0; i < base_size; i++) {
	    search_stack[i].cmp_to_best_path = 0;
	    search_stack[i].in_best_path = true;
	    PathInfo path_info;
	    path_info.splitting_element = search_stack[i].split_element;
	    path_info.certificate_index = search_stack[i].certificate_index;
	    path_info.eqref_hash = search_stack[i].eqref_hash;
	    path_info.subcertificate_length = search_stack[i].subcertificate_length;
	    best_path_info.push_back(path_info);
	    //fprintf(stdout, "%u ", search_stack[i].split_element);
	  }
	  certificate_best_path = certificate_current_path;
	  //fprintf(stdout, "\n"); fflush(stdout);
	  /*
	   * Backtrack to the previous level
	   */
	  p.level--;
	  continue;
	}

      {
	//fprintf(stderr, "BAUT ");
	/*
	 * Equal to the previous best path
	 */
#if defined(BLISS_CONSISTENCY_CHECKS)
	/* Verify that the automorphism is correctly built */
	for(unsigned int i = 0; i < N; i++)
	  assert(best_path_automorphism[i] ==
		 p.elements[best_path_labeling[i]]);
#endif
	
#if defined(BLISS_VERIFY_AUTOMORPHISMS)
	/* Verify that it really is an automorphism */
	assert(is_automorphism(best_path_automorphism));
#endif
      
	unsigned int gca_level_with_first = 0;
	for(unsigned int i = search_stack.size(); i > 0; i--) {
	  if((int)first_path_info[gca_level_with_first].splitting_element !=
	     search_stack[gca_level_with_first].split_element)
	    break;
	  gca_level_with_first++;
	}
	assert(gca_level_with_first < p.level);

	unsigned int gca_level_with_best = 0;
	for(unsigned int i = search_stack.size(); i > 0; i--) {
	  if((int)best_path_info[gca_level_with_best].splitting_element !=
	     search_stack[gca_level_with_best].split_element)
	    break;
	  gca_level_with_best++;
	}
	assert(gca_level_with_best < p.level);

	if(opt_use_long_prune)
	  {
	    long_prune_add_automorphism(best_path_automorphism);
	  }
	    
	/*
	 * Update orbit information
	 */
	update_orbit_information(best_path_orbits, best_path_automorphism);

	/*
	 * Update orbit information
	 */
	const unsigned int nof_old_orbits = first_path_orbits.nof_orbits();
	update_orbit_information(first_path_orbits, best_path_automorphism);
	if(nof_old_orbits != first_path_orbits.nof_orbits())
	  {
	    /* Some orbits were merged */
	    /* Report automorphism */
	    if(report_hook)
	      (*report_hook)(report_user_param,
			     get_nof_vertices(),
			     best_path_automorphism);
	    /* Update statistics */
	    stats.nof_generators++;
	  }
	  
	/*
	 * Compute backjumping level
	 */
	unsigned int backjumping_level = p.level - 1;
	if(!first_path_orbits.is_minimal_representative(search_stack[gca_level_with_first].split_element))
	  {
	    backjumping_level = gca_level_with_first;
	    /*fprintf(stderr, "bj1: %u %u\n", p.level, backjumping_level);*/
	  }
	else
	  {
	    assert(!best_path_orbits.is_minimal_representative(search_stack[gca_level_with_best].split_element));
	    backjumping_level = gca_level_with_best;
	    /*fprintf(stderr, "bj2: %u %u\n", p.level, backjumping_level);*/
	  }
	/* Backtrack */
	search_stack.resize(backjumping_level + 1);
	p.level = backjumping_level;
	continue;
      }
    }


  if(opt_use_long_prune)
    {
      long_prune_deallocate();
    }


}




void AbstractGraph::find_automorphisms(Stats &stats,
				       void (*hook)(void *user_param,
						    unsigned int n,
						    const unsigned int *aut),
				       void *user_param)
{
  report_hook = hook;
  report_user_param = user_param;

  search(false, stats);

  if(first_path_labeling)
    {
      free(first_path_labeling);
      first_path_labeling = 0;
    }
  if(best_path_labeling)
    {
      free(best_path_labeling);
      best_path_labeling = 0;
    }
}


const unsigned int *
AbstractGraph::canonical_form(Stats &stats,
			      void (*hook)(void *user_param,
					   unsigned int n,
					   const unsigned int *aut),
			      void *user_param)
{

  report_hook = hook;
  report_user_param = user_param;

  search(true, stats);

  return best_path_labeling;
}




/*-------------------------------------------------------------------------
 *
 * Routines for directed graphs
 *
 *-------------------------------------------------------------------------*/

Digraph::Vertex::Vertex()
{
  color = 0;
}


Digraph::Vertex::~Vertex()
{
  ;
}


void Digraph::Vertex::add_edge_to(const unsigned int other_vertex)
{
  edges_out.push_back(other_vertex);
}


void Digraph::Vertex::add_edge_from(const unsigned int other_vertex)
{
  edges_in.push_back(other_vertex);
}


void Digraph::Vertex::remove_duplicate_edges(bool * const duplicate_array)
{
  for(std::vector<unsigned int>::iterator iter = edges_out.begin();
      iter != edges_out.end(); )
    {
      const unsigned int dest_vertex = *iter;
      if(duplicate_array[dest_vertex] == true)
	{
	  /* A duplicate edge found! */
	  iter = edges_out.erase(iter);
	}
      else
	{
	  /* Not seen earlier, mark as seen */
	  duplicate_array[dest_vertex] = true;
	  iter++;
	}
    }

  /* Clear duplicate_array */
  for(std::vector<unsigned int>::iterator iter = edges_out.begin();
      iter != edges_out.end();
      iter++)
    {
      duplicate_array[*iter] = false;
    }

  for(std::vector<unsigned int>::iterator iter = edges_in.begin();
      iter != edges_in.end(); )
    {
      const unsigned int dest_vertex = *iter;
      if(duplicate_array[dest_vertex] == true)
	{
	  /* A duplicate edge found! */
	  iter = edges_in.erase(iter);
	}
      else
	{
	  /* Not seen earlier, mark as seen */
	  duplicate_array[dest_vertex] = true;
	  iter++;
	}
    }

  /* Clear duplicate_array */
  for(std::vector<unsigned int>::iterator iter = edges_in.begin();
      iter != edges_in.end();
      iter++)
    {
      duplicate_array[*iter] = false;
    }
}


/**
 * Sort the edges entering and leaving the vertex according to
 * the vertex number of the other edge end.
 * Time complexity: O(e log(e)), where e is the number of edges
 * entering/leaving the vertex.
 */
void Digraph::Vertex::sort_edges()
{
  std::sort(edges_in.begin(), edges_in.end());
  std::sort(edges_out.begin(), edges_out.end());
}





/*-------------------------------------------------------------------------
 *
 * Constructor and destructor for directed graphs
 *
 *-------------------------------------------------------------------------*/


Digraph::Digraph(const unsigned int nof_vertices)
{
  vertices.resize(nof_vertices);
}


Digraph::~Digraph()
{
  ;
}


unsigned int Digraph::add_vertex(const unsigned int color)
{
  const unsigned int new_vertex_num = vertices.size();
  vertices.resize(new_vertex_num + 1);
  vertices.back().color = color;
  return new_vertex_num;
}


void Digraph::add_edge(const unsigned int vertex1, const unsigned int vertex2)
{
  assert(vertex1 < vertices.size());
  assert(vertex2 < vertices.size());
  vertices[vertex1].add_edge_to(vertex2);
  vertices[vertex2].add_edge_from(vertex1);
}


void Digraph::change_color(const unsigned int vertex,
			   const unsigned int new_color)
{
  assert(vertex < vertices.size());
  vertices[vertex].color = new_color;
}


void Digraph::sort_edges()
{
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    vertices[i].sort_edges();
}





Digraph *Digraph::permute(const unsigned int *perm) const
{
  Digraph *g = new Digraph(get_nof_vertices());
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      const Vertex &v = vertices[i];
      g->change_color(perm[i], v.color);
      for(std::vector<unsigned int>::const_iterator ei = v.edges_out.begin();
	  ei != v.edges_out.end();
	  ei++)
	{
	  g->add_edge(perm[i], perm[*ei]);
	}
    }
  g->sort_edges();
  return g;
}





/*-------------------------------------------------------------------------
 *
 * Print graph in graphviz format
 *
 *-------------------------------------------------------------------------*/


void Digraph::write_dot(const char * const file_name)
{
  FILE *fp = fopen(file_name, "w");
  if(fp)
    write_dot(fp);
  fclose(fp);
}


void Digraph::write_dot(FILE * const fp)
{
  remove_duplicate_edges();

  fprintf(fp, "digraph g {\n");

  unsigned int vnum = 0;
  for(std::vector<Vertex>::iterator vi = vertices.begin();
      vi != vertices.end();
      vi++, vnum++)
    {
      Vertex &v = *vi;
      fprintf(fp, "v%u [label=\"%u:%u\"];\n", vnum, vnum, v.color);
      for(std::vector<unsigned int>::const_iterator ei = v.edges_out.begin();
	  ei != v.edges_out.end();
	  ei++)
	{
	  const unsigned int vnum2 = *ei;
	  fprintf(fp, "v%u -> v%u\n", vnum, vnum2);
	}
    }

  fprintf(fp, "}\n");
}


void Digraph::remove_duplicate_edges()
{
  bool *duplicate_array = (bool*)calloc(vertices.size(), sizeof(bool));

  for(std::vector<Vertex>::iterator vi = vertices.begin();
      vi != vertices.end();
      vi++)
    {
#if defined(BLISS_EXPENSIVE_CONSISTENCY_CHECKS)
      for(unsigned int i = 0; i < vertices.size(); i++)
	assert(duplicate_array[i] == false);
#endif
      Vertex &v = *vi;
      v.remove_duplicate_edges(duplicate_array);
    }

  free(duplicate_array);
}





/*-------------------------------------------------------------------------
 *
 * Get a hash value for the graph.
 *
 *-------------------------------------------------------------------------*/
unsigned int Digraph::get_hash()
{
  remove_duplicate_edges();
  sort_edges();

  UintSeqHash h;

  h.update(get_nof_vertices());

  /* Hash the color of each vertex */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      h.update(vertices[i].color);
    }

  /* Hash the edges */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges_out.begin();
	  ei != v.edges_out.end();
	  ei++)
	{
	  h.update(i);
	  h.update(*ei);
	}
    }

  return h.get_value();
}



/*-------------------------------------------------------------------------
 *
 * Read directed graph in the DIMACS format.
 * Returns 0 if an error occurred.
 *
 *-------------------------------------------------------------------------*/

Digraph *Digraph::read_dimacs(FILE * const fp, FILE * const errstr)
{
  Digraph *g = 0;
  unsigned int nof_vertices;
  unsigned int nof_edges;
  unsigned int line_num = 1;
  int c;

  const bool verbose = false;
  FILE * const verbstr = stdout;
  
  /* Read comments and the problem definition line */
  while(1)
    {
      c = getc(fp);
      if(c == 'c')
	{
	  /* A comment, ignore the rest of the line */
	  while((c = getc(fp)) != '\n')
	    {
	      if(c == EOF) {
		if(errstr)
		  fprintf(errstr, "error in line %u: not in DIMACS format\n",
			  line_num);
		goto error_exit;
	      }
	    }
	  line_num++;
	  continue;
	}
      if(c == 'p')
	{
	  /* The problem definition line */
	  if(fscanf(fp, " edge %u %u\n", &nof_vertices, &nof_edges) != 2)
	    {
	      if(errstr)
		fprintf(errstr, "error in line %u: not in DIMACS format\n",
			line_num);
	      goto error_exit;
	    }
	  line_num++;
	  break;
	}
      if(errstr)
	fprintf(errstr, "error in line %u: not in DIMACS format\n", line_num);
      goto error_exit;
    }
  
  if(nof_vertices <= 0)
    {
      if(errstr)
	fprintf(errstr, "error: no vertices\n");
      goto error_exit;
    }
  if(verbose)
    {
      fprintf(verbstr, "Instance has %d vertices and %d edges\n",
	      nof_vertices, nof_edges);
      fflush(verbstr);
    }

  g = new Digraph(nof_vertices);

  //
  // Read vertex colors
  //
  if(verbose)
    {
      fprintf(verbstr, "Reading vertex colors...\n");
      fflush(verbstr);
    }
  while(1)
    {
      c = getc(fp);
      if(c != 'n')
	{
	  ungetc(c, fp);
	  break;
	}
      ungetc(c, fp);
      unsigned int vertex;
      unsigned int color;
      if(fscanf(fp, "n %u %u\n", &vertex, &color) != 2)
	{
	  if(errstr)
	    fprintf(errstr, "error in line %u: not in DIMACS format\n",
		    line_num);
	  goto error_exit;
	}
      if(!((vertex >= 1) && (vertex <= nof_vertices)))
	{
	  if(errstr)
	    fprintf(errstr,
		    "error in line %u: vertex %u not in range [1,...%u]\n",
		    line_num, vertex, nof_vertices);
	  goto error_exit;
	}
      line_num++;
      g->change_color(vertex - 1, color);
    }
  if(verbose)
    {
      fprintf(verbstr, "Done\n");
      fflush(verbstr);
    }

  //
  // Read edges
  //
  if(verbose)
    {
      fprintf(verbstr, "Reading edges...\n");
      fflush(verbstr);
    }
  for(unsigned i = 0; i < nof_edges; i++)
    {
      unsigned int from, to;
      if(fscanf(fp, "e %u %u\n", &from, &to) != 2)
	{
	  if(errstr)
	    fprintf(errstr, "error in line %u: not in DIMACS format\n",
		    line_num);
	  goto error_exit;
	}
      if(!((from >= 1) && (from <= nof_vertices)))
	{
	  if(errstr)
	    fprintf(errstr,
		    "error in line %u: vertex %u not in range [1,...%u]\n",
		    line_num, from, nof_vertices);
	  goto error_exit;
	}
      if(!((to >= 1) && (to <= nof_vertices)))
	{
	  if(errstr)
	    fprintf(errstr,
		    "error in line %u: vertex %u not in range [1,...%u]\n",
		    line_num, to, nof_vertices);
	  goto error_exit;
	}
      line_num++;
      g->add_edge(from-1, to-1);
    }
  if(verbose)
    {
      fprintf(verbstr, "Done\n");
      fflush(verbstr);
    }
  
  return g;

 error_exit:
  if(g)
    delete g;
  return 0;
}





void Digraph::write_dimacs(FILE * const fp)
{
  remove_duplicate_edges();
  sort_edges();

  /* First count the total number of edges */
  unsigned int nof_edges = 0;
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      nof_edges += vertices[i].edges_out.size();
    }

  /* Output the "header" line */
  fprintf(fp, "p edge %u %u\n", get_nof_vertices(), nof_edges);

  /* Print the color of each vertex */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      fprintf(fp, "n %u %u\n", i+1, v.color);
      /*
      if(v.color != 0)
	{
	  fprintf(fp, "n %u %u\n", i+1, v.color);
	}
      */
    }

  /* Print the edges */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges_out.begin();
	  ei != v.edges_out.end();
	  ei++)
	{
	  fprintf(fp, "e %u %u\n", i+1, (*ei)+1);
	}
    }
}








/*-------------------------------------------------------------------------
 *
 * Partition independent invariants
 *
 *-------------------------------------------------------------------------*/


unsigned int Digraph::vertex_color_invariant(Digraph * const g,
					     const unsigned int vertex)
{
  BLISS_ASSERT(vertex < g->vertices.size());
  return g->vertices[vertex].color;
}


unsigned int Digraph::indegree_invariant(Digraph * const g,
					 const unsigned int vertex)
{
  BLISS_ASSERT(vertex < g->vertices.size());
  return g->vertices[vertex].get_nof_edges_in();
}


unsigned int Digraph::outdegree_invariant(Digraph * const g,
					  const unsigned int vertex)
{
  BLISS_ASSERT(vertex < g->vertices.size());
  return g->vertices[vertex].get_nof_edges_out();
}



unsigned int Digraph::selfloop_invariant(Digraph * const g,
					 const unsigned int v)
{
  /* VERY INEFFICIENT! */
  BLISS_ASSERT(v < g->vertices.size());
  Vertex &vertex = g->vertices[v];
  for(std::vector<unsigned int>::iterator ei = vertex.edges_out.begin();
      ei != vertex.edges_out.end();
      ei++)
    {
      if(*ei == v)
	return 1;
    }
  return 0;
}





/*-------------------------------------------------------------------------
 *
 * Refine the partition p according to a partition independent invariant
 *
 *-------------------------------------------------------------------------*/

bool Digraph::refine_according_to_invariant(unsigned int (*inv)(Digraph * const g, unsigned int v))
{
  bool refined = false;

  for(Partition::Cell *cell = p.first_cell; cell; )
    {
      assert(cell->max_ival == 0);
      assert(cell->max_ival_count == 0);
      
      Partition::Cell * const next_cell = cell->next;
      
      if(cell->is_unit())
	{
	  cell = next_cell;
	  continue;
	}

      const unsigned int *ep = p.elements + cell->first;
      for(unsigned int i = cell->length; i > 0; i--, ep++)
	{
	  unsigned int ival = inv(this, *ep);
	  p.invariant_values[*ep] = ival;
	  if(ival > cell->max_ival) {
	    cell->max_ival = ival;
	    cell->max_ival_count = 1;
	  }
	  else if(ival == cell->max_ival) {
	    cell->max_ival_count++;
	  }
	}
      Partition::Cell * const last_new_cell = p.zplit_cell(cell, true);
      refined = (last_new_cell != cell);
      cell = next_cell;
    }

  return refined;
}





/*-------------------------------------------------------------------------
 *
 * Split the neighbourhood of a cell according to the equitable invariant
 *
 *-------------------------------------------------------------------------*/


void Digraph::split_neighbourhood_of_cell(Partition::Cell * const cell)
{
  BLISS_ASSERT(neighbour_heap.is_empty());
  BLISS_ASSERT(cell->length > 1);
  
  if(compute_eqref_hash)
    {
      eqref_hash.update(cell->first);
      eqref_hash.update(cell->length);
    }

  const unsigned int *ep = p.elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--)
    {
      const Vertex &v = vertices[*ep++];
      
      std::vector<unsigned int>::const_iterator ei = v.edges_out.begin();
      for(unsigned int j = v.get_nof_edges_out(); j > 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->is_unit())
	    continue;
	  const unsigned int ival = ++p.invariant_values[dest_vertex];
	  if(ival > neighbour_cell->max_ival) {
	    neighbour_cell->max_ival = ival;
	    neighbour_cell->max_ival_count = 1;
	    if(ival == 1)
	      neighbour_heap.insert(neighbour_cell->first);
	  }
	  else if(ival == neighbour_cell->max_ival) {
	    neighbour_cell->max_ival_count++;
	  }
	}
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell * const neighbour_cell =
	p.element_to_cell_map[p.elements[start]];
      BLISS_ASSERT(neighbour_cell->first == start);
      BLISS_ASSERT(neighbour_cell->length > 1);
      BLISS_ASSERT(neighbour_cell->max_ival >= 1);
      BLISS_ASSERT(neighbour_cell->max_ival_count >= 1);
      
      if(compute_eqref_hash)
	{
	  eqref_hash.update(neighbour_cell->first);
	  eqref_hash.update(neighbour_cell->length);
	  eqref_hash.update(neighbour_cell->max_ival);
	  eqref_hash.update(neighbour_cell->max_ival_count);
	}
      Partition::Cell * const last_new_cell = p.zplit_cell(neighbour_cell, true);
      /* Update hash */
      const Partition::Cell *c = neighbour_cell;
      while(1)
	{
	  if(compute_eqref_hash)
	    {
	      eqref_hash.update(c->first);
	      eqref_hash.update(c->length);
	    }
	  if(c == last_new_cell)
	    break;
	  c = c->next;
	}
    }

  if(cell->in_splitting_queue)
    {
      return;
    }

  ep = p.elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--)
    {
      const Vertex &v = vertices[*ep++];

      std::vector<unsigned int>::const_iterator ei = v.edges_in.begin();
      for(unsigned int j = v.get_nof_edges_in(); j > 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Partition::Cell * const neighbour_cell =
	    p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->is_unit())
	    continue;
	  const unsigned int ival = ++p.invariant_values[dest_vertex];
	  if(ival > neighbour_cell->max_ival)
	    {
	      neighbour_cell->max_ival = ival;
	      neighbour_cell->max_ival_count = 1;
	      if(ival == 1)
		neighbour_heap.insert(neighbour_cell->first);
	    }
	  else if(ival == neighbour_cell->max_ival) {
	    neighbour_cell->max_ival_count++;
	  }
	}
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell * const neighbour_cell =
	p.element_to_cell_map[p.elements[start]];
      BLISS_ASSERT(neighbour_cell->first == start);
      BLISS_ASSERT(neighbour_cell->length > 1);
      BLISS_ASSERT(neighbour_cell->max_ival >= 1);
      BLISS_ASSERT(neighbour_cell->max_ival_count >= 1);
      
      if(compute_eqref_hash)
	{
	  eqref_hash.update(neighbour_cell->first);
	  eqref_hash.update(neighbour_cell->length);
	  eqref_hash.update(neighbour_cell->max_ival);
	  eqref_hash.update(neighbour_cell->max_ival_count);
	}

      Partition::Cell * const last_new_cell = p.zplit_cell(neighbour_cell, true);
      /* Update hash */
      const Partition::Cell *c = neighbour_cell;
      while(1)
	{
	  if(compute_eqref_hash)
	    {
	      eqref_hash.update(c->first);
	      eqref_hash.update(c->length);
	    }
	  if(c == last_new_cell)
	    break;
	  c = c->next;
	}
    }
}


bool Digraph::split_neighbourhood_of_unit_cell(Partition::Cell * const unit_cell)
{
  BLISS_ASSERT(neighbour_heap.is_empty());

  BLISS_ASSERT(unit_cell->is_unit());
  BLISS_ASSERT(p.element_to_cell_map[p.elements[unit_cell->first]] ==
	       unit_cell);
  BLISS_ASSERT(p.in_pos[p.elements[unit_cell->first]] ==
               p.elements + unit_cell->first);

  if(compute_eqref_hash)
    {
      eqref_hash.update(0x87654321);
      eqref_hash.update(unit_cell->first);
      eqref_hash.update(1);
    }

  const Vertex &v = vertices[p.elements[unit_cell->first]];

  /*
   * Phase 1
   * Refine neighbours according to the edges that leave the vertex v
   */
  std::vector<unsigned int>::const_iterator ei = v.edges_out.begin();
  for(unsigned int j = v.get_nof_edges_out(); j > 0; j--)
    {
      const unsigned int dest_vertex = *ei++;
      Partition::Cell * const neighbour_cell =
	p.element_to_cell_map[dest_vertex];
      BLISS_ASSERT(*p.in_pos[dest_vertex] == dest_vertex);
   
      if(neighbour_cell->is_unit()) {
	if(in_search) {
	  neighbour_heap.insert(neighbour_cell->first);
	}
	continue;
      }
      if(neighbour_cell->max_ival_count == 0)
	{
	  neighbour_heap.insert(neighbour_cell->first);
	}
      neighbour_cell->max_ival_count++;
      BLISS_ASSERT(neighbour_cell->max_ival_count <= neighbour_cell->length);
      
      unsigned int * const swap_position =
	p.elements + neighbour_cell->first + neighbour_cell->length -
	neighbour_cell->max_ival_count;
      BLISS_ASSERT(p.in_pos[dest_vertex] <= swap_position);
      *p.in_pos[dest_vertex] = *swap_position;
      p.in_pos[*swap_position] = p.in_pos[dest_vertex];
      *swap_position = dest_vertex;
      p.in_pos[dest_vertex] = swap_position;
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell *neighbour_cell =
	p.element_to_cell_map[p.elements[start]];

#if defined(BLISS_CONSISTENCY_CHECKS)
      assert(neighbour_cell->first == start);
      if(neighbour_cell->is_unit()) {
	assert(neighbour_cell->max_ival_count == 0);
      } else {
	assert(neighbour_cell->max_ival_count > 0);
	assert(neighbour_cell->max_ival_count <= neighbour_cell->length);
      }
#endif

      if(compute_eqref_hash)
	{
	  eqref_hash.update(neighbour_cell->first);
	  eqref_hash.update(neighbour_cell->length);
	  eqref_hash.update(neighbour_cell->max_ival_count);
	}

      if(neighbour_cell->length > 1 &&
	 neighbour_cell->max_ival_count != neighbour_cell->length) {
	
	p.consistency_check();

	Partition::Cell * const new_cell =
	  p.aux_split_in_two(neighbour_cell,
			     neighbour_cell->length -
			     neighbour_cell->max_ival_count);
	unsigned int *ep = p.elements + new_cell->first;
	unsigned int * const lp = p.elements+new_cell->first+new_cell->length;
	while(ep < lp) {
	  BLISS_ASSERT(p.in_pos[*ep] == ep);
	  p.element_to_cell_map[*ep] = new_cell;
	  ep++;
	}
	neighbour_cell->max_ival_count = 0;

	p.consistency_check();

	if(compute_eqref_hash)
	  {
	    eqref_hash.update(neighbour_cell->first);
	    eqref_hash.update(neighbour_cell->length);
	    eqref_hash.update(0);
	    eqref_hash.update(new_cell->first);
	    eqref_hash.update(new_cell->length);
	    eqref_hash.update(1);
	  }

	/* Add cells in splitting_queue */
	BLISS_ASSERT(!new_cell->in_splitting_queue);
	if(neighbour_cell->in_splitting_queue) {
	  /* Both cells must be included in splitting_queue in order
	     to have refinement to equitable partition */
	  p.add_in_splitting_queue(new_cell);
	} else {
	  Partition::Cell *min_cell, *max_cell;
	  if(neighbour_cell->length <= new_cell->length) {
	    min_cell = neighbour_cell;
	    max_cell = new_cell;
	  } else {
	    min_cell = new_cell;
	    max_cell = neighbour_cell;
	  }
	  /* Put the smaller cell in splitting_queue */
	   p.add_in_splitting_queue(min_cell);
	  if(max_cell->is_unit()) {
	    /* Put the "larger" cell also in splitting_queue */
	    p.add_in_splitting_queue(max_cell);
	  }
	}
	/* Update pointer for certificate generation */
	neighbour_cell = new_cell;
      }
      else
	{
	  neighbour_cell->max_ival_count = 0;
	}
      
      /*
       * Build certificate if required
       */
      if(in_search)
	{
	  for(unsigned int i = neighbour_cell->first,
		j = neighbour_cell->length;
	      j > 0;
	      j--, i++)
	    {
	      if(refine_compare_certificate)
		{
		  if(refine_equal_to_first)
		    {
		      if(refine_current_path_certificate_index >=
			 refine_first_path_subcertificate_end)
			{
			  refine_equal_to_first = false;
			}
		      else if(certificate_first_path[refine_current_path_certificate_index] != unit_cell->first)
			{
			  refine_equal_to_first = false;
			}
		      else if(certificate_first_path[refine_current_path_certificate_index+1] != i)
			{
			  refine_equal_to_first = false;
			}
		    }
		  if(refine_cmp_to_best == 0)
		    {
		      if(refine_current_path_certificate_index >=
			 refine_best_path_subcertificate_end)
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first > certificate_best_path[refine_current_path_certificate_index])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first < certificate_best_path[refine_current_path_certificate_index])
			{
			  refine_cmp_to_best = -1;
			}
		      else if(i > certificate_best_path[refine_current_path_certificate_index+1])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(i < certificate_best_path[refine_current_path_certificate_index+1])
			{
			  refine_cmp_to_best = -1;
			}
		    }
		  if((refine_equal_to_first == false) &&
		     (refine_cmp_to_best < 0))
		    goto worse_exit;
		}
	      certificate_current_path[refine_current_path_certificate_index++] = unit_cell->first;
	      certificate_current_path[refine_current_path_certificate_index++] = i;
	    }
	} /* if(in_search) */
    } /* while(!neighbour_heap.is_empty()) */

  /*
   * Phase 2
   * Refine neighbours according to the edges that enter the vertex v
   */
  ei = v.edges_in.begin();
  for(unsigned int j = v.get_nof_edges_in(); j > 0; j--)
    {
      const unsigned int dest_vertex = *ei++;
      Partition::Cell * const neighbour_cell =
	p.element_to_cell_map[dest_vertex];
      BLISS_ASSERT(*p.in_pos[dest_vertex] == dest_vertex);
      
      if(neighbour_cell->is_unit()) {
	if(in_search) {
	  neighbour_heap.insert(neighbour_cell->first);
	}
	continue;
      }
      if(neighbour_cell->max_ival_count == 0)
	{
	  neighbour_heap.insert(neighbour_cell->first);
	}
      neighbour_cell->max_ival_count++;
      BLISS_ASSERT(neighbour_cell->max_ival_count <= neighbour_cell->length);

      unsigned int * const swap_position =
	p.elements + neighbour_cell->first + neighbour_cell->length -
	neighbour_cell->max_ival_count;
      BLISS_ASSERT(p.in_pos[dest_vertex] <= swap_position);
      *p.in_pos[dest_vertex] = *swap_position;
      p.in_pos[*swap_position] = p.in_pos[dest_vertex];
      *swap_position = dest_vertex;
      p.in_pos[dest_vertex] = swap_position;
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell *neighbour_cell =
	p.element_to_cell_map[p.elements[start]];

#if defined(BLISS_CONSISTENCY_CHECKS)
      assert(neighbour_cell->first == start);
      if(neighbour_cell->is_unit()) {
	assert(neighbour_cell->max_ival_count == 0);
      } else {
	assert(neighbour_cell->max_ival_count > 0);
	assert(neighbour_cell->max_ival_count <= neighbour_cell->length);
      }
#endif

      if(compute_eqref_hash)
	{
	  eqref_hash.update(neighbour_cell->first);
	  eqref_hash.update(neighbour_cell->length);
	  eqref_hash.update(neighbour_cell->max_ival_count);
	}

      if(neighbour_cell->length > 1 &&
	 neighbour_cell->max_ival_count != neighbour_cell->length) {

	p.consistency_check();

	Partition::Cell * const new_cell =
	  p.aux_split_in_two(neighbour_cell,
			     neighbour_cell->length -
			     neighbour_cell->max_ival_count);
	unsigned int *ep = p.elements + new_cell->first;
	unsigned int * const lp = p.elements+new_cell->first+new_cell->length;
	while(ep < lp) {
	  BLISS_ASSERT(p.in_pos[*ep] == ep);
	  p.element_to_cell_map[*ep] = new_cell;
	  ep++;
	}
	neighbour_cell->max_ival_count = 0;

	p.consistency_check();

	if(compute_eqref_hash)
	  {
	    eqref_hash.update(neighbour_cell->first);
	    eqref_hash.update(neighbour_cell->length);
	    eqref_hash.update(0);
	    eqref_hash.update(new_cell->first);
	    eqref_hash.update(new_cell->length);
	    eqref_hash.update(1);
	  }

	/* Add cells in splitting_queue */
	BLISS_ASSERT(!new_cell->in_splitting_queue);
	if(neighbour_cell->in_splitting_queue) {
	  /* Both cells must be included in splitting_queue in order
	     to have refinement to equitable partition */
	  p.add_in_splitting_queue(new_cell);
	} else {
	  Partition::Cell *min_cell, *max_cell;
	  if(neighbour_cell->length <= new_cell->length) {
	    min_cell = neighbour_cell;
	    max_cell = new_cell;
	  } else {
	    min_cell = new_cell;
	    max_cell = neighbour_cell;
	  }
	  /* Put the smaller cell in splitting_queue */
	   p.add_in_splitting_queue(min_cell);
	  if(max_cell->is_unit()) {
	    /* Put the "larger" cell also in splitting_queue */
	    p.add_in_splitting_queue(max_cell);
	  }
	}
	/* Update pointer for certificate generation */
	neighbour_cell = new_cell;
      }
      else
	{
	  neighbour_cell->max_ival_count = 0;
	}
      
      /*
       * Build certificate if required
       */
      if(in_search)
	{
	  for(unsigned int i = neighbour_cell->first,
		j = neighbour_cell->length;
	      j > 0;
	      j--, i++)
	    {
	      if(refine_compare_certificate)
		{
		  if(refine_equal_to_first)
		    {
		      if(refine_current_path_certificate_index >=
			 refine_first_path_subcertificate_end)
			{
			  refine_equal_to_first = false;
			}
		      else if(certificate_first_path[refine_current_path_certificate_index] != i)
			{
			  refine_equal_to_first = false;
			}
		      else if(certificate_first_path[refine_current_path_certificate_index+1] != unit_cell->first)
			{
			  refine_equal_to_first = false;
			}
		    }
		  if(refine_cmp_to_best == 0)
		    {
		      if(refine_current_path_certificate_index >=
			 refine_best_path_subcertificate_end)
			{
			  refine_cmp_to_best = 1;
			}
		      else if(i > certificate_best_path[refine_current_path_certificate_index])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(i < certificate_best_path[refine_current_path_certificate_index])
			{
			  refine_cmp_to_best = -1;
			}
		      else if(unit_cell->first > certificate_best_path[refine_current_path_certificate_index+1])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first < certificate_best_path[refine_current_path_certificate_index+1])
			{
			  refine_cmp_to_best = -1;
			}
		    }
		  if((refine_equal_to_first == false) &&
		     (refine_cmp_to_best < 0))
		    goto worse_exit;
		}
	      certificate_current_path[refine_current_path_certificate_index++] = i;
	      certificate_current_path[refine_current_path_certificate_index++] = unit_cell->first;
	    }
	} /* if(in_search) */
    } /* while(!neighbour_heap.is_empty()) */

  return false;

 worse_exit:
  /* Clear neighbour heap */
  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell * const neighbour_cell = p.element_to_cell_map[p.elements[start]];
      neighbour_cell->max_ival_count = 0;
    }
  return true;
}





/*-------------------------------------------------------------------------
 *
 * Check whether the current partition p is equitable.
 * Performance: very slow, use only for debugging purposes.
 *
 *-------------------------------------------------------------------------*/

bool Digraph::is_equitable() const
{
  const unsigned int N = get_nof_vertices();
  if(N == 0)
    return true;

  std::vector<unsigned int> first_count = std::vector<unsigned int>(N, 0);
  std::vector<unsigned int> other_count = std::vector<unsigned int>(N, 0);

  /*
   * Check equitabledness w.r.t. outgoing edges
   */
  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->is_unit())
	continue;

      unsigned int *ep = p.elements + cell->first;
      const Vertex &first_vertex = vertices[*ep++];

      /* Count outgoing edges of the first vertex for cells */
      for(std::vector<unsigned int>::const_iterator ei =
	    first_vertex.edges_out.begin();
	  ei != first_vertex.edges_out.end();
	  ei++)
	{
	  first_count[p.element_to_cell_map[*ei]->first]++;
	}

      /* Count and compare outgoing edges of the other vertices */
      for(unsigned int i = cell->length; i > 1; i--)
	{
	  const Vertex &vertex = vertices[*ep++];
	  for(std::vector<unsigned int>::const_iterator ei =
		vertex.edges_out.begin();
	      ei != vertex.edges_out.end();
	      ei++)
	    {
	      other_count[p.element_to_cell_map[*ei]->first]++;
	    }
	  for(Partition::Cell *cell2 = p.first_cell;
	      cell2;
	      cell2 = cell2->next)
	    {
	      if(first_count[cell2->first] != other_count[cell2->first])
		{
		  /* Not equitable */
		  return false;
		}
	      other_count[cell2->first] = 0;
	    }
	}
      /* Reset first_count */
      for(unsigned int i = 0; i < N; i++)
	first_count[i] = 0;
    }


  /*
   * Check equitabledness w.r.t. incoming edges
   */
  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->is_unit())
	continue;

      unsigned int *ep = p.elements + cell->first;
      const Vertex &first_vertex = vertices[*ep++];

      /* Count incoming edges of the first vertex for cells */
      for(std::vector<unsigned int>::const_iterator ei =
	    first_vertex.edges_in.begin();
	  ei != first_vertex.edges_in.end();
	  ei++)
	{
	  first_count[p.element_to_cell_map[*ei]->first]++;
	}

      /* Count and compare incoming edges of the other vertices */
      for(unsigned int i = cell->length; i > 1; i--)
	{
	  const Vertex &vertex = vertices[*ep++];
	  for(std::vector<unsigned int>::const_iterator ei =
		vertex.edges_in.begin();
	      ei != vertex.edges_in.end();
	      ei++)
	    {
	      other_count[p.element_to_cell_map[*ei]->first]++;
	    }
	  for(Partition::Cell *cell2 = p.first_cell;
	      cell2;
	      cell2 = cell2->next)
	    {
	      if(first_count[cell2->first] != other_count[cell2->first])
		{
		  /* Not equitable */
		  return false;
		}
	      other_count[cell2->first] = 0;
	    }
	}
      /* Reset first_count */
      for(unsigned int i = 0; i < N; i++)
	first_count[i] = 0;
    }
  return true;
}





/*-------------------------------------------------------------------------
 *
 * Build the initial equitable partition
 *
 *-------------------------------------------------------------------------*/

void Digraph::make_initial_equitable_partition()
{
  refine_according_to_invariant(&vertex_color_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_according_to_invariant(&outdegree_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_according_to_invariant(&indegree_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_according_to_invariant(&selfloop_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_to_equitable();
  //p.print_signature(stderr); fprintf(stderr, "\n");
}





/*-------------------------------------------------------------------------
 *
 * Find the next cell to be splitted
 *
 *-------------------------------------------------------------------------*/

Partition::Cell *Digraph::find_next_cell_to_be_splitted(Partition::Cell *cell)
{
  assert(!p.is_discrete());
  switch(sh) {
  case shs_f:
    return sh_first(cell);
  case shs_fs:
    return sh_first_smallest(cell);
  case shs_fl:
    return sh_first_largest(cell);
  case shs_fm:
    return sh_first_max_neighbours(cell);
  case shs_fsm:
    return sh_first_smallest_max_neighbours(cell);
  case shs_flm:
    return sh_first_largest_max_neighbours(cell);
  default:
    assert(false && "Unknown splitting heuristics");
  }
}

/** \internal
 * A splitting heuristic.
 * Returns the first nonsingleton cell in the current partition.
 * The argument \a cell is ignored.
 */
Partition::Cell *Digraph::sh_first(Partition::Cell *cell)
{
  return p.first_nonsingleton_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first smallest nonsingleton cell in the current partition.
 * The argument \a cell is ignored.
 */
Partition::Cell *Digraph::sh_first_smallest(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  unsigned int best_size = UINT_MAX;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
      if(cell->length < best_size)
	{
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first largest nonsingleton cell in the current partition.
 * The argument \a cell is ignored.
 */
Partition::Cell *Digraph::sh_first_largest(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  unsigned int best_size = 0;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
      if(cell->length > best_size)
	{
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Assumes that the max_ival fields of the cells are all 0.
 */
Partition::Cell *Digraph::sh_first_max_neighbours(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  int best_value = -1;
  KStack<Partition::Cell*> neighbour_cells_visited;
  neighbour_cells_visited.init(get_nof_vertices());
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
      int value = 0;
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei;

      ei = v.edges_in.begin();
      for(unsigned int j = v.get_nof_edges_in(); j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}

      ei = v.edges_out.begin();
      for(unsigned int j = v.get_nof_edges_out(); j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      
      if(value > best_value)
	{
	  best_value = value;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first smallest nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Assumes that the max_ival fields of the cells are all 0.
 */
Partition::Cell *Digraph::sh_first_smallest_max_neighbours(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  int best_value = -1;
  unsigned int best_size = UINT_MAX;
  KStack<Partition::Cell*> neighbour_cells_visited;
  neighbour_cells_visited.init(get_nof_vertices());
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
	
      int value = 0;
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei;

      ei = v.edges_in.begin();
      for(unsigned int j = v.get_nof_edges_in(); j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}

      ei = v.edges_out.begin();
      for(unsigned int j = v.get_nof_edges_out(); j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}

      if((value > best_value) ||
	 (value == best_value && cell->length < best_size))
	{
	  best_value = value;
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first largest nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Assumes that the max_ival fields of the cells are all 0.
 */
Partition::Cell *Digraph::sh_first_largest_max_neighbours(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  int best_value = -1;
  unsigned int best_size = 0;
  KStack<Partition::Cell*> neighbour_cells_visited;
  neighbour_cells_visited.init(get_nof_vertices());
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
	
      int value = 0;
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei;

      ei = v.edges_in.begin();
      for(unsigned int j = v.get_nof_edges_in(); j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}

      ei = v.edges_out.begin();
      for(unsigned int j = v.get_nof_edges_out(); j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}

      if((value > best_value) ||
	 (value == best_value && cell->length > best_size))
	{
	  best_value = value;
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}




/*------------------------------------------------------------------------
 *
 * Initialize the certificate size and memory
 *
 *-------------------------------------------------------------------------*/

void Digraph::initialize_certificate()
{
  certificate_size = 0;
  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->length > 1) {
	certificate_size +=
	  vertices[p.elements[cell->first]].get_nof_edges_in() * 2 *
	  cell->length;
	certificate_size +=
	  vertices[p.elements[cell->first]].get_nof_edges_out() * 2 *
	  cell->length;
      }
    }
  certificate_index = 0;

  certificate_current_path.clear();
  certificate_current_path.resize(certificate_size);
  certificate_first_path.clear();
  certificate_best_path.clear();
}



bool Digraph::is_automorphism(unsigned int * const perm)
{
  std::set<unsigned int, std::less<unsigned int> > edges1;
  std::set<unsigned int, std::less<unsigned int> > edges2;

  bool result = true;

  for(unsigned int i = 0; i < vertices.size(); i++)
    {
      Vertex &v1 = vertices[i];
      Vertex &v2 = vertices[perm[i]];

      edges1.clear();
      for(std::vector<unsigned int>::iterator ei = v1.edges_in.begin();
	  ei != v1.edges_in.end();
	  ei++)
	edges1.insert(perm[*ei]);
      edges2.clear();
      for(std::vector<unsigned int>::iterator ei = v2.edges_in.begin();
	  ei != v2.edges_in.end();
	  ei++)
	edges2.insert(*ei);
      if(!(edges1 == edges2))
	{
	  result = false;
	  goto done;
	}

      edges1.clear();
      for(std::vector<unsigned int>::iterator ei = v1.edges_out.begin();
	  ei != v1.edges_out.end();
	  ei++)
	edges1.insert(perm[*ei]);
      edges2.clear();
      for(std::vector<unsigned int>::iterator ei = v2.edges_out.begin();
	  ei != v2.edges_out.end();
	  ei++)
	edges2.insert(*ei);
      if(!(edges1 == edges2))
	{
	  result = false;
	  goto done;
	}
    }

 done:

  return result;
}





/*-------------------------------------------------------------------------
 *
 * Routines for undirected graphs
 *
 *-------------------------------------------------------------------------*/

Graph::Vertex::Vertex()
{
  color = 0;
  nof_edges = 0;
}


Graph::Vertex::~Vertex()
{
  ;
}


void Graph::Vertex::add_edge(const unsigned int other_vertex)
{
  edges.push_back(other_vertex);
  nof_edges++;
  BLISS_ASSERT(nof_edges == edges.size());
}


void Graph::Vertex::remove_duplicate_edges(bool * const duplicate_array)
{
  for(std::vector<unsigned int>::iterator iter = edges.begin();
      iter != edges.end(); )
    {
      const unsigned int dest_vertex = *iter;
      if(duplicate_array[dest_vertex] == true)
	{
	  /* A duplicate edge found! */
	  iter = edges.erase(iter);
	  nof_edges--;
	  BLISS_ASSERT(nof_edges == edges.size());
	}
      else
	{
	  /* Not seen earlier, mark as seen */
	  duplicate_array[dest_vertex] = true;
	  iter++;
	}
    }

  /* Clear duplicate_array */
  for(std::vector<unsigned int>::iterator iter = edges.begin();
      iter != edges.end();
      iter++)
    {
      duplicate_array[*iter] = false;
    }
}


/**
 * Sort the edges leaving the vertex according to
 * the vertex number of the other edge end.
 * Time complexity: O(e log(e)), where e is the number of edges
 * leaving the vertex.
 */
void Graph::Vertex::sort_edges()
{
  std::sort(edges.begin(), edges.end());
}





/*-------------------------------------------------------------------------
 *
 * Constructor and destructor for undirected graphs
 *
 *-------------------------------------------------------------------------*/


Graph::Graph(const unsigned int nof_vertices)
{
  vertices.resize(nof_vertices);
  sh = shs_flm;
}


Graph::~Graph()
{
  ;
}


unsigned int Graph::add_vertex(const unsigned int color)
{
  const unsigned int vertex_num = vertices.size();
  vertices.resize(vertex_num + 1);
  vertices.back().color = color;
  return vertex_num;
}


void Graph::add_edge(const unsigned int vertex1, const unsigned int vertex2)
{
  //fprintf(stderr, "(%u,%u) ", vertex1, vertex2);
  assert(vertex1 < vertices.size());
  assert(vertex2 < vertices.size());
  vertices[vertex1].add_edge(vertex2);
  vertices[vertex2].add_edge(vertex1);
}


void Graph::change_color(const unsigned int vertex,
			 const unsigned int color)
{
  assert(vertex < vertices.size());
  vertices[vertex].color = color;
}





/*-------------------------------------------------------------------------
 *
 * Read graph in the DIMACS format.
 * Returns 0 if an error occurred.
 *
 *-------------------------------------------------------------------------*/

Graph *Graph::read_dimacs(FILE * const fp, FILE * const errstr)
{
  Graph *g = 0;
  unsigned int nof_vertices;
  unsigned int nof_edges;
  unsigned int line_num = 1;
  int c;

  const bool verbose = false;
  FILE * const verbstr = stdout;
  
  /* Read comments and the problem definition line */
  while(1)
    {
      c = getc(fp);
      if(c == 'c')
	{
	  /* A comment, ignore the rest of the line */
	  while((c = getc(fp)) != '\n')
	    {
	      if(c == EOF)
		{
		  if(errstr)
		    fprintf(errstr,
			    "error in line %u: not in DIMACS format\n",
			    line_num);
		  goto error_exit;
		}
	    }
	  line_num++;
	  continue;
	}
      if(c == 'p')
	{
	  /* The problem definition line */
	  if(fscanf(fp, " edge %u %u\n", &nof_vertices, &nof_edges) != 2)
	    {
	      if(errstr)
		fprintf(errstr, "error in line %u: not in DIMACS format\n",
			line_num);
	      goto error_exit;
	    }
	  line_num++;
	  break;
	}
      if(errstr)
	fprintf(errstr, "error in line %u: not in DIMACS format\n", line_num);
      goto error_exit;
    }
  
  if(nof_vertices <= 0)
    {
      if(errstr)
	fprintf(errstr, "error: no vertices\n");
      goto error_exit;
    }
  if(verbose)
    {
      fprintf(verbstr, "Instance has %d vertices and %d edges\n",
	      nof_vertices, nof_edges);
      fflush(verbstr);
    }

  g = new Graph(nof_vertices);

  //
  // Read vertex colors
  //
  if(verbose)
    {
      fprintf(verbstr, "Reading vertex colors...\n");
      fflush(verbstr);
    }
  while(1)
    {
      c = getc(fp);
      if(c != 'n')
	{
	  ungetc(c, fp);
	  break;
	}
      ungetc(c, fp);
      unsigned int vertex;
      unsigned int color;
      if(fscanf(fp, "n %u %u\n", &vertex, &color) != 2)
	{
	  if(errstr)
	    fprintf(errstr, "error in line %u: not in DIMACS format\n",
		    line_num);
	  goto error_exit;
	}
      if(!((vertex >= 1) && (vertex <= nof_vertices)))
	{
	  if(errstr)
	    fprintf(errstr,
		    "error in line %u: vertex %u not in range [1,...,%u]\n",
		    line_num, vertex, nof_vertices);
	  goto error_exit;
	}
      line_num++;
      g->change_color(vertex - 1, color);
    }
  if(verbose)
    {
      fprintf(verbstr, "Done\n");
      fflush(verbstr);
    }

  //
  // Read edges
  //
  if(verbose)
    {
      fprintf(verbstr, "Reading edges...\n");
      fflush(verbstr);
    }
  for(unsigned i = 0; i < nof_edges; i++)
    {
      unsigned int from, to;
      if(fscanf(fp, "e %u %u\n", &from, &to) != 2)
	{
	  if(errstr)
	    fprintf(errstr, "error in line %u: not in DIMACS format\n",
		    line_num);
	  goto error_exit;
	}
      if(!((from >= 1) && (from <= nof_vertices)))
	{
	  if(errstr)
	    fprintf(errstr,
		    "error in line %u: vertex %u not in range [1,...,%u]\n",
		    line_num, from, nof_vertices);
	  goto error_exit;
	}
      if(!((to >= 1) && (to <= nof_vertices)))
	{
	  if(errstr)
	    fprintf(errstr,
		    "error in line %u: vertex %u not in range [1,...,%u]\n",
		    line_num, to, nof_vertices);
	  goto error_exit;
	}
      line_num++;
      g->add_edge(from-1, to-1);
    }
  if(verbose)
    {
      fprintf(verbstr, "Done\n");
      fflush(verbstr);
    }

  return g;

 error_exit:
  if(g)
    delete g;
  return 0;

}


void Graph::write_dimacs(FILE * const fp)
{
  remove_duplicate_edges();
  sort_edges();

  /* First count the total number of edges */
  unsigned int nof_edges = 0;
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_i = *ei;
	  if(dest_i < i)
	    continue;
	  nof_edges++;
	}
    }

  /* Output the "header" line */
  fprintf(fp, "p edge %u %u\n", get_nof_vertices(), nof_edges);

  /* Print the color of each vertex */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      fprintf(fp, "n %u %u\n", i+1, v.color);
      /*
      if(v.color != 0)
	{
	  fprintf(fp, "n %u %u\n", i+1, v.color);
	}
      */
    }

  /* Print the edges */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_i = *ei;
	  if(dest_i < i)
	    continue;
	  fprintf(fp, "e %u %u\n", i+1, dest_i+1);
	}
    }
}



void Graph::sort_edges()
{
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    vertices[i].sort_edges();
}


int Graph::cmp(Graph *g)
{
  /* Compare the numbers of vertices */
  if(get_nof_vertices() < g->get_nof_vertices())
    return -1;
  if(get_nof_vertices() > g->get_nof_vertices())
    return 1;
  /* Compare vertex colors and degrees */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      if(vertices[i].color < g->vertices[i].color)
	return -1;
      if(vertices[i].color > g->vertices[i].color)
	return 1;
      if(vertices[i].nof_edges < g->vertices[i].nof_edges)
	return -1;
      if(vertices[i].nof_edges > g->vertices[i].nof_edges)
	return 1;
    }
  /* Compare edges */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      Vertex &gv = g->vertices[i];
      v.sort_edges();
      gv.sort_edges();
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      std::vector<unsigned int>::const_iterator gei = gv.edges.begin();
      while(ei != v.edges.end())
	{
	  if(*ei < *gei)
	    return -1;
	  if(*ei > *gei)
	    return 1;
	  ei++;
	  gei++;
	}
    }
  return 0;
}


Graph *Graph::permute(const unsigned int *perm) const
{
  Graph *g = new Graph(get_nof_vertices());
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      const Vertex &v = vertices[i];
      Vertex &permuted_v = g->vertices[perm[i]];
      permuted_v.color = v.color;
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_v = *ei;
	  permuted_v.add_edge(perm[dest_v]);
	}
      permuted_v.sort_edges();
    }
  return g;
}





/*-------------------------------------------------------------------------
 *
 * Print graph in graphviz format
 *
 *-------------------------------------------------------------------------*/


void Graph::write_dot(const char * const file_name)
{
  FILE *fp = fopen(file_name, "w");
  if(fp)
    write_dot(fp);
  fclose(fp);
}


void Graph::write_dot(FILE * const fp)
{
  remove_duplicate_edges();

  fprintf(fp, "graph g {\n");

  unsigned int vnum = 0;
  for(std::vector<Vertex>::iterator vi = vertices.begin();
      vi != vertices.end();
      vi++, vnum++)
    {
      Vertex &v = *vi;
      fprintf(fp, "v%u [label=\"%u:%u\"];\n", vnum, vnum, v.color);
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int vnum2 = *ei;
	  if(vnum2 > vnum)
	    fprintf(fp, "v%u -- v%u\n", vnum, vnum2);
	}
    }

  fprintf(fp, "}\n");
}





/*-------------------------------------------------------------------------
 *
 * Get a hash value for the graph.
 *
 *-------------------------------------------------------------------------*/

unsigned int Graph::get_hash()
{
  remove_duplicate_edges();
  sort_edges();

  UintSeqHash h;

  h.update(get_nof_vertices());

  /* Hash the color of each vertex */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      h.update(vertices[i].color);
    }

  /* Hash the edges */
  for(unsigned int i = 0; i < get_nof_vertices(); i++)
    {
      Vertex &v = vertices[i];
      for(std::vector<unsigned int>::const_iterator ei = v.edges.begin();
	  ei != v.edges.end();
	  ei++)
	{
	  const unsigned int dest_i = *ei;
	  if(dest_i < i)
	    continue;
	  h.update(i);
	  h.update(dest_i);
	}
    }

  return h.get_value();
}





void Graph::remove_duplicate_edges()
{
  bool *duplicate_array = (bool*)calloc(vertices.size(), sizeof(bool));

  for(std::vector<Vertex>::iterator vi = vertices.begin();
      vi != vertices.end();
      vi++)
    {
#if defined(BLISS_EXPENSIVE_CONSISTENCY_CHECKS)
      for(unsigned int i = 0; i < vertices.size(); i++)
	assert(duplicate_array[i] == false);
#endif
      Vertex &v = *vi;
      v.remove_duplicate_edges(duplicate_array);
    }

  free(duplicate_array);
}





/*-------------------------------------------------------------------------
 *
 * Partition independent invariants
 *
 *-------------------------------------------------------------------------*/


/*
 * Return the color of the vertex.
 * Time complexity: O(1)
 */
unsigned int Graph::vertex_color_invariant(Graph * const g,
					   const unsigned int v)
{
  BLISS_ASSERT(v < g->vertices.size());
  return g->vertices[v].color;
}


/*
 * Return the degree of the vertex.
 * Time complexity: O(1)
 */
unsigned int Graph::degree_invariant(Graph * const g,
				     const unsigned int v)
{
  BLISS_ASSERT(v < g->vertices.size());
  BLISS_ASSERT(g->vertices[v].edges.size() ==
	       g->vertices[v].nof_edges);
  return g->vertices[v].nof_edges;
}


/*
 * Return 1 if the vertex v has a self-loop, 0 otherwise
 * Time complexity: O(E_v), where E_v is the number of edges leaving v
 */
unsigned int Graph::selfloop_invariant(Graph *g, unsigned int v)
{
  BLISS_ASSERT(v < g->vertices.size());
  BLISS_ASSERT(g->vertices[v].edges.size() ==
	       g->vertices[v].nof_edges);
  Vertex &vertex = g->vertices[v];
  for(std::vector<unsigned int>::const_iterator ei = vertex.edges.begin();
      ei != vertex.edges.end();
      ei++)
    {
      if(*ei == v)
	return 1;
    }
  return 0;
}






/*-------------------------------------------------------------------------
 *
 * Refine the partition p according to a partition independent invariant
 *
 *-------------------------------------------------------------------------*/

bool Graph::refine_according_to_invariant(unsigned int (*inv)(Graph * const g,
							      unsigned int v))
{
  bool refined = false;

  for(Partition::Cell *cell = p.first_cell; cell; )
    {
      BLISS_ASSERT(cell->max_ival == 0);
      BLISS_ASSERT(cell->max_ival_count == 0);
      
      Partition::Cell * const next_cell = cell->next;

      if(cell->is_unit())
	{
	  cell = next_cell;
	  continue;
	}
      
      const unsigned int *ep = p.elements + cell->first;
      for(unsigned int i = cell->length; i > 0; i--, ep++)
	{
	  unsigned int ival = inv(this, *ep);
	  p.invariant_values[*ep] = ival;
	  if(ival > cell->max_ival)
	    {
	      cell->max_ival = ival;
	      cell->max_ival_count = 1;
	    }
	  else if(ival == cell->max_ival)
	    {
	      cell->max_ival_count++;
	    }
	}
      Partition::Cell * const last_new_cell = p.zplit_cell(cell, true);
      refined = (last_new_cell != cell);
      cell = next_cell;
    }

  return refined;
}










/*-------------------------------------------------------------------------
 *
 * Split the neighbourhood of a cell according to the equitable invariant
 *
 *-------------------------------------------------------------------------*/


void Graph::split_neighbourhood_of_cell(Partition::Cell * const cell)
{
  BLISS_ASSERT(neighbour_heap.is_empty());
  BLISS_ASSERT(cell->length > 1);

  if(compute_eqref_hash)
    {
      eqref_hash.update(cell->first);
      eqref_hash.update(cell->length);
    }

  unsigned int *ep = p.elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--)
    {
      const Vertex &v = vertices[*ep++];
      
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      for(unsigned int j = v.nof_edges; j != 0; j--)
	{
	  const unsigned int dest_vertex = *ei++;
	  Partition::Cell * const neighbour_cell =
	    p.element_to_cell_map[dest_vertex];
	  if(neighbour_cell->is_unit())
	    continue;
	  const unsigned int ival = ++p.invariant_values[dest_vertex];
	  if(ival > neighbour_cell->max_ival)
	    {
	      neighbour_cell->max_ival = ival;
	      neighbour_cell->max_ival_count = 1;
	      if(ival == 1)
		neighbour_heap.insert(neighbour_cell->first);
	    }
	  else if(ival == neighbour_cell->max_ival) {
	    neighbour_cell->max_ival_count++;
	  }
	}
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell * const neighbour_cell = p.element_to_cell_map[p.elements[start]];
      BLISS_ASSERT(neighbour_cell->first == start);
      BLISS_ASSERT(neighbour_cell->length > 1);
      BLISS_ASSERT(neighbour_cell->max_ival >= 1);
      BLISS_ASSERT(neighbour_cell->max_ival_count >= 1);
      
      if(compute_eqref_hash)
	{
	  eqref_hash.update(neighbour_cell->first);
	  eqref_hash.update(neighbour_cell->length);
	  eqref_hash.update(neighbour_cell->max_ival);
	  eqref_hash.update(neighbour_cell->max_ival_count);
	}

      Partition::Cell * const last_new_cell = p.zplit_cell(neighbour_cell, true);
      /* Update hash */
      const Partition::Cell *c = neighbour_cell;
      while(1)
	{
	  if(compute_eqref_hash)
	    {
	      eqref_hash.update(c->first);
	      eqref_hash.update(c->length);
	    }
	  if(c == last_new_cell)
	    break;
	  c = c->next;
	}
    }
}


bool Graph::split_neighbourhood_of_unit_cell(Partition::Cell * const unit_cell)
{
  BLISS_ASSERT(neighbour_heap.is_empty());

  BLISS_ASSERT(unit_cell->is_unit());
  BLISS_ASSERT(p.element_to_cell_map[p.elements[unit_cell->first]] ==
	       unit_cell);
  BLISS_ASSERT(p.in_pos[p.elements[unit_cell->first]] ==
               p.elements + unit_cell->first);

  if(compute_eqref_hash)
    {
      eqref_hash.update(0x87654321);
      eqref_hash.update(unit_cell->first);
      eqref_hash.update(1);
    }

  const Vertex &v = vertices[p.elements[unit_cell->first]];

  std::vector<unsigned int>::const_iterator ei = v.edges.begin();
  for(unsigned int j = v.nof_edges; j > 0; j--)
    {
      const unsigned int dest_vertex = *ei++;
      Partition::Cell * const neighbour_cell =
	p.element_to_cell_map[dest_vertex];
      BLISS_ASSERT(*p.in_pos[dest_vertex] == dest_vertex);
      
      if(neighbour_cell->is_unit()) {
	if(in_search) {
	  neighbour_heap.insert(neighbour_cell->first);
	}
	continue;
      }
      if(neighbour_cell->max_ival_count == 0)
	{
	  neighbour_heap.insert(neighbour_cell->first);
	}
      neighbour_cell->max_ival_count++;
      BLISS_ASSERT(neighbour_cell->max_ival_count <= neighbour_cell->length);
      unsigned int * const swap_position =
	p.elements + neighbour_cell->first + neighbour_cell->length -
	neighbour_cell->max_ival_count;
      BLISS_ASSERT(p.in_pos[dest_vertex] <= swap_position);
      *p.in_pos[dest_vertex] = *swap_position;
      p.in_pos[*swap_position] = p.in_pos[dest_vertex];
      *swap_position = dest_vertex;
      p.in_pos[dest_vertex] = swap_position;
    }

  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell *neighbour_cell =
	p.element_to_cell_map[p.elements[start]];

#if defined(BLISS_CONSISTENCY_CHECKS)
      assert(neighbour_cell->first == start);
      if(neighbour_cell->is_unit()) {
	assert(neighbour_cell->max_ival_count == 0);
      } else {
	assert(neighbour_cell->max_ival_count > 0);
	assert(neighbour_cell->max_ival_count <= neighbour_cell->length);
      }
#endif

      if(compute_eqref_hash)
	{
	  eqref_hash.update(neighbour_cell->first);
	  eqref_hash.update(neighbour_cell->length);
	  eqref_hash.update(neighbour_cell->max_ival_count);
	}

      if(neighbour_cell->length > 1 &&
	 neighbour_cell->max_ival_count != neighbour_cell->length) {

	p.consistency_check();

	Partition::Cell * const new_cell =
	  p.aux_split_in_two(neighbour_cell,
			     neighbour_cell->length -
			     neighbour_cell->max_ival_count);
	unsigned int *ep = p.elements + new_cell->first;
	unsigned int * const lp = p.elements+new_cell->first+new_cell->length;
	while(ep < lp)
	  {
	    BLISS_ASSERT(p.in_pos[*ep] == ep);
	    p.element_to_cell_map[*ep] = new_cell;
	    ep++;
	  }
	neighbour_cell->max_ival_count = 0;

	p.consistency_check();

	/* update hash */
	if(compute_eqref_hash)
	  {
	    eqref_hash.update(neighbour_cell->first);
	    eqref_hash.update(neighbour_cell->length);
	    eqref_hash.update(0);
	    eqref_hash.update(new_cell->first);
	    eqref_hash.update(new_cell->length);
	    eqref_hash.update(1);
	  }

	/* Add cells in splitting_queue */
	BLISS_ASSERT(!new_cell->in_splitting_queue);
	if(neighbour_cell->in_splitting_queue) {
	  /* Both cells must be included in splitting_queue in order
	     to have refinement to equitable partition */
	  p.add_in_splitting_queue(new_cell);
	} else {
	  Partition::Cell *min_cell, *max_cell;
	  if(neighbour_cell->length <= new_cell->length) {
	    min_cell = neighbour_cell;
	    max_cell = new_cell;
	  } else {
	    min_cell = new_cell;
	    max_cell = neighbour_cell;
	  }
	  /* Put the smaller cell in splitting_queue */
	  p.add_in_splitting_queue(min_cell);
	  if(max_cell->is_unit()) {
	    /* Put the "larger" cell also in splitting_queue */
	    p.add_in_splitting_queue(max_cell);
	  }
	}
	/* Update pointer for certificate generation */
	neighbour_cell = new_cell;
      }
      else
	{
	  neighbour_cell->max_ival_count = 0;
	}

      /*
       * Build certificate if required
       */
      if(in_search)
	{
	  for(unsigned int i = neighbour_cell->first,
		j = neighbour_cell->length;
	      j > 0;
	      j--, i++)
	    {
	      if(refine_compare_certificate)
		{
		  if(refine_equal_to_first)
		    {
		      if(refine_current_path_certificate_index >=
			 refine_first_path_subcertificate_end)
			{
			  refine_equal_to_first = false;
			}
		      else if(certificate_first_path[refine_current_path_certificate_index] != unit_cell->first)
			{
			  refine_equal_to_first = false;
			}
		      else if(certificate_first_path[refine_current_path_certificate_index+1] != i)
			{
			  refine_equal_to_first = false;
			}
		    }
		  if(refine_cmp_to_best == 0)
		    {
		      if(refine_current_path_certificate_index >=
			 refine_best_path_subcertificate_end)
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first > certificate_best_path[refine_current_path_certificate_index])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(unit_cell->first < certificate_best_path[refine_current_path_certificate_index])
			{
			  refine_cmp_to_best = -1;
			}
		      else if(i > certificate_best_path[refine_current_path_certificate_index+1])
			{
			  refine_cmp_to_best = 1;
			}
		      else if(i < certificate_best_path[refine_current_path_certificate_index+1])
			{
			  refine_cmp_to_best = -1;
			}
		    }
		  if((refine_equal_to_first == false) &&
		     (refine_cmp_to_best < 0))
		    goto worse_exit;
		}
	      certificate_current_path[refine_current_path_certificate_index++] = unit_cell->first;
	      certificate_current_path[refine_current_path_certificate_index++] = i;
	    }
	} /* if(in_search) */
    } /* while(!neighbour_heap.is_empty()) */
  
  return false;

 worse_exit:
  /* Clear neighbour heap */
  while(!neighbour_heap.is_empty())
    {
      const unsigned int start = neighbour_heap.remove();
      Partition::Cell * const neighbour_cell =
	p.element_to_cell_map[p.elements[start]];
      neighbour_cell->max_ival_count = 0;
    }
  return true;
}





/*-------------------------------------------------------------------------
 *
 * Check whether the current partition p is equitable.
 * Performance: very slow, use only for debugging purposes.
 *
 *-------------------------------------------------------------------------*/

bool Graph::is_equitable() const
{
  const unsigned int N = get_nof_vertices();
  if(N == 0)
    return true;

  std::vector<unsigned int> first_count = std::vector<unsigned int>(N, 0);
  std::vector<unsigned int> other_count = std::vector<unsigned int>(N, 0);

  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->is_unit())
	continue;
      
      unsigned int *ep = p.elements + cell->first;
      const Vertex &first_vertex = vertices[*ep++];

      /* Count how many edges lead from the first vertex to
       * the neighbouring cells */
      for(std::vector<unsigned int>::const_iterator ei =
	    first_vertex.edges.begin();
	  ei != first_vertex.edges.end();
	  ei++)
	{
	  first_count[p.element_to_cell_map[*ei]->first]++;
	}

      /* Count and compare to the edges of the other vertices */
      for(unsigned int i = cell->length; i > 1; i--)
	{
	  const Vertex &vertex = vertices[*ep++];
	  for(std::vector<unsigned int>::const_iterator ei =
		vertex.edges.begin();
	      ei != vertex.edges.end();
	      ei++)
	    {
	      other_count[p.element_to_cell_map[*ei]->first]++;
	    }
	  for(Partition::Cell *cell2 = p.first_cell;
	      cell2;
	      cell2 = cell2->next)
	    {
	      if(first_count[cell2->first] != other_count[cell2->first])
		{
		  /* Not equitable */
		  return false;
		}
	      other_count[cell2->first] = 0;
	    }
	}
      /* Reset first_count */
      for(unsigned int i = 0; i < N; i++)
	first_count[i] = 0;
    }
  return true;
}





/*-------------------------------------------------------------------------
 *
 * Build the initial equitable partition
 *
 *-------------------------------------------------------------------------*/

void Graph::make_initial_equitable_partition()
{
  refine_according_to_invariant(&vertex_color_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_according_to_invariant(&degree_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_according_to_invariant(&selfloop_invariant);
  p.clear_splitting_queue();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_to_equitable();
  //p.print_signature(stderr); fprintf(stderr, "\n");

  refine_to_equitable();
  //p.print_signature(stderr); fprintf(stderr, "\n");

}





/*-------------------------------------------------------------------------
 *
 * Find the next cell to be splitted
 *
 *-------------------------------------------------------------------------*/


Partition::Cell *Graph::find_next_cell_to_be_splitted(Partition::Cell *cell)
{
  assert(!p.is_discrete());
  switch(sh) {
  case shs_f:
    return sh_first(cell);
  case shs_fs:
    return sh_first_smallest(cell);
  case shs_fl:
    return sh_first_largest(cell);
  case shs_fm:
    return sh_first_max_neighbours(cell);
  case shs_fsm:
    return sh_first_smallest_max_neighbours(cell);
  case shs_flm:
    return sh_first_largest_max_neighbours(cell);
  default:
    assert(false && "Unknown splitting heuristics");
  }
}

/** \internal
 * A splitting heuristic.
 * Returns the first nonsingleton cell in the current partition.
 * The argument \A cell is ignored.
 */
Partition::Cell *Graph::sh_first(Partition::Cell *cell)
{
  return p.first_nonsingleton_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first smallest nonsingleton cell in the current partition.
 * The argument \a cell is ignored.
 */
Partition::Cell *Graph::sh_first_smallest(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  unsigned int best_size = UINT_MAX;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
      if(cell->length < best_size)
	{
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first largest nonsingleton cell in the current partition.
 * The argument \a cell is ignored.
 */
Partition::Cell *Graph::sh_first_largest(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  unsigned int best_size = 0;
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
      if(cell->length > best_size)
	{
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Assumes that the max_ival fields of the cells are all 0.
 */
Partition::Cell *Graph::sh_first_max_neighbours(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  int best_value = -1;
  KStack<Partition::Cell*> neighbour_cells_visited;
  neighbour_cells_visited.init(get_nof_vertices());
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
	
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      int value = 0;
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      if(value > best_value)
	{
	  best_value = value;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first smallest nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Assumes that the max_ival fields of the cells are all 0.
 */
Partition::Cell *Graph::sh_first_smallest_max_neighbours(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  int best_value = -1;
  unsigned int best_size = UINT_MAX;
  KStack<Partition::Cell*> neighbour_cells_visited;
  neighbour_cells_visited.init(get_nof_vertices());
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
	
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      int value = 0;
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      if((value > best_value) ||
	 (value == best_value && cell->length < best_size))
	{
	  best_value = value;
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}

/** \internal
 * A splitting heuristic.
 * Returns the first largest nonsingleton cell with max number of neighbouring
 * nonsingleton cells.
 * Assumes that the partition p is equitable.
 * Assumes that the max_ival fields of the cells are all 0.
 */
Partition::Cell *Graph::sh_first_largest_max_neighbours(Partition::Cell *cell)
{
  Partition::Cell *best_cell = 0;
  int best_value = -1;
  unsigned int best_size = 0;
  KStack<Partition::Cell*> neighbour_cells_visited;
  neighbour_cells_visited.init(get_nof_vertices());
  for(cell = p.first_nonsingleton_cell; cell; cell = cell->next_nonsingleton)
    {
      BLISS_ASSERT(cell->length > 1);
	
      const Vertex &v = vertices[p.elements[cell->first]];
      std::vector<unsigned int>::const_iterator ei = v.edges.begin();
      for(unsigned int j = v.nof_edges; j > 0; j--)
	{
	  Partition::Cell * const neighbour_cell = p.element_to_cell_map[*ei++];
	  if(neighbour_cell->is_unit())
	    continue;
	  neighbour_cell->max_ival++;
	  if(neighbour_cell->max_ival == 1)
	    neighbour_cells_visited.push(neighbour_cell);
	}
      int value = 0;
      while(!neighbour_cells_visited.is_empty())
	{
	  Partition::Cell * const neighbour_cell = neighbour_cells_visited.pop();
	  if(neighbour_cell->max_ival != neighbour_cell->length)
	    value++;
	  neighbour_cell->max_ival = 0;
	}
      if((value > best_value) ||
	 (value == best_value && cell->length > best_size))
	{
	  best_value = value;
	  best_size = cell->length;
	  best_cell = cell;
	}
    }
  BLISS_ASSERT(best_cell != 0);
  return best_cell;
}





/*-------------------------------------------------------------------------
 *
 * Initialize the certificate size and memory
 *
 *-------------------------------------------------------------------------*/

void Graph::initialize_certificate()
{
  certificate_size = 0;
  for(Partition::Cell *cell = p.first_cell; cell; cell = cell->next)
    {
      if(cell->length > 1) {
	certificate_size +=
	  vertices[p.elements[cell->first]].nof_edges * 2 * cell->length;
      }
    }
  certificate_index = 0;

  certificate_current_path.clear();
  certificate_current_path.resize(certificate_size);
  certificate_first_path.clear();
  certificate_best_path.clear();
}





/*-------------------------------------------------------------------------
 *
 * Check whether perm is an automorphism
 *
 *-------------------------------------------------------------------------*/

bool Graph::is_automorphism(unsigned int * const perm)
{
  std::set<unsigned int, std::less<unsigned int> > edges1;
  std::set<unsigned int, std::less<unsigned int> > edges2;

  bool result = true;

  for(unsigned int i = 0; i < vertices.size(); i++)
    {
      Vertex &v1 = vertices[i];
      edges1.clear();
      for(std::vector<unsigned int>::iterator ei = v1.edges.begin();
	  ei != v1.edges.end();
	  ei++)
	edges1.insert(perm[*ei]);
      
      Vertex &v2 = vertices[perm[i]];
      edges2.clear();
      for(std::vector<unsigned int>::iterator ei = v2.edges.begin();
	  ei != v2.edges.end();
	  ei++)
	edges2.insert(*ei);

      if(!(edges1 == edges2))
	{
	  result = false;
	  goto done;
	}
    }

 done:

  return result;
}







}
