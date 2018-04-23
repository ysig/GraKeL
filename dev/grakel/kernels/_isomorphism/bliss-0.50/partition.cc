#include <assert.h>
#include <vector>
#include <list>
#include "graph.hh"
#include "partition.hh"

/*
 * Copyright (c) Tommi Junttila
 * Released under the GNU General Public License version 2.
 */


namespace bliss {

static const bool should_not_happen = false;



Partition::Partition()
{
  N = 0;
  elements = 0;
  in_pos = 0;
  invariant_values = 0;
  cells = 0;
  free_cells = 0;
  element_to_cell_map = 0;
  graph = 0;
  /* Initialize a distribution count sorting array. */
  for(unsigned int i = 0; i < 256; i++)
    dcs_count[i] = 0;
}



Partition::~Partition()
{
  if(elements) {
    free(elements); elements = 0; }
  if(cells) {
    free(cells); cells = 0; }
  if(element_to_cell_map) {
    free(element_to_cell_map); element_to_cell_map = 0; }
  if(in_pos) {
    free(in_pos); in_pos = 0; }
  if(invariant_values) {
    free(invariant_values); invariant_values = 0; }
  N = 0;
}



void Partition::init(const unsigned int M)
{
  assert(M > 0);
  N = M;

  if(elements)
    free(elements);
  elements = (unsigned int*)malloc(N * sizeof(unsigned int));
  for(unsigned int i = 0; i < N; i++)
    elements[i] = i;

  if(in_pos)
    free(in_pos);
  in_pos = (unsigned int**)malloc(N * sizeof(unsigned int*));
  for(unsigned int i = 0; i < N; i++)
    in_pos[i] = elements + i;

  if(invariant_values)
    free(invariant_values);
  invariant_values = (unsigned int*)malloc(N * sizeof(unsigned int));
  for(unsigned int i = 0; i < N; i++)
    invariant_values[i] = 0;

  if(cells)
    free(cells);
  cells = (Cell*)malloc(N * sizeof(Cell));

  cells[0].first = 0;
  cells[0].length = N;
  cells[0].max_ival = 0;
  cells[0].max_ival_count = 0;
  cells[0].in_splitting_queue = false;
  cells[0].in_neighbour_heap = false;
  cells[0].next = 0;
  cells[0].prev_next_ptr = &first_cell;
  cells[0].next_nonsingleton = 0;
  cells[0].prev_nonsingleton = 0;
  cells[0].split_level = 0;
  first_cell = &cells[0];
  if(N == 1)
    first_nonsingleton_cell = 0;
  else
    first_nonsingleton_cell = &cells[0];

  for(unsigned int i = 1; i < N; i++)
    {
      cells[i].first = 0;
      cells[i].length = 0;
      cells[i].max_ival = 0;
      cells[i].max_ival_count = 0;
      cells[i].in_splitting_queue = false;
      cells[i].in_neighbour_heap = false;
      cells[i].next = (i < N-1)?&cells[i+1]:0;
      cells[i].prev_next_ptr = (i == 1)?&free_cells:&(cells[i-1].next);
      cells[i].next_nonsingleton = 0;
      cells[i].prev_nonsingleton = 0;
    }
  if(N > 1)
    free_cells = &cells[1];
  else
    free_cells = 0;

  if(element_to_cell_map)
    free(element_to_cell_map);
  element_to_cell_map = (Cell **)malloc(N * sizeof(Cell *));
  for(unsigned int i = 0; i < N; i++)
    element_to_cell_map[i] = first_cell;

  splitting_queue.init(N);
  refinement_stack.init(N);
  level = 0;

  /* Reset the main backtracking stack */
  bt_stack.clear();
}

/*
 * For debugging purposes only.
 */
void Partition::consistency_check()
{
#if defined(BLISS_CONSISTENCY_CHECKS)
  for(unsigned int i = 0; i < N; i++)
    {
      assert(*in_pos[i] == i);
      assert(in_pos[elements[i]] == elements + i);
      assert((unsigned int)(in_pos[i] - elements) >=
	     element_to_cell_map[i]->first);
      assert((unsigned int)(in_pos[i] - elements) <
	     element_to_cell_map[i]->first + element_to_cell_map[i]->length);
    }
  for(const Cell *cell = first_cell; cell; cell = cell->next)
    {
      assert(cell->prev_next_ptr && *(cell->prev_next_ptr) == cell);
    }
  const bool do_print = false;
  if(do_print)
    {
      fprintf(stderr, "\nRef stack: ");
      for(unsigned int j = 0; j < refinement_stack.size(); j++)
	{
	  const RefInfo i = refinement_stack.element_at(j);
	  fprintf(stderr, "f%u,%d,%d ",
	      i.split_cell_first,
		  i.prev_nonsingleton_first,
		  i.next_nonsingleton_first);
	}
      fprintf(stderr, "\n");
      for(const Cell *cell = first_nonsingleton_cell; cell;
	  cell = cell->next_nonsingleton)
	{
	  fprintf(stderr, "%u:%u->", cell->first, cell->length);
	  if(cell->next_nonsingleton)
	    assert(cell->first < cell->next_nonsingleton->first);
	}
      fprintf(stderr, "\n");
    }
  const Cell *next_nonsingleton = first_nonsingleton_cell;
  const Cell *prev_nonsingleton = 0;
  if(next_nonsingleton)
    assert(next_nonsingleton->prev_nonsingleton == 0);
  for(const Cell *cell = first_cell; cell; cell = cell->next)
    {
      assert(!cell->next || cell->next->prev_next_ptr == &(cell->next));
      if(cell->length > 1)
	{
	  if(do_print)
	    fprintf(stderr, "%u:%u=>", cell->first, cell->length);
	  assert(cell == next_nonsingleton); 
	  assert(cell->prev_nonsingleton == prev_nonsingleton);
	  next_nonsingleton = cell->next_nonsingleton;
	  prev_nonsingleton = cell;
	  if(next_nonsingleton)
	    assert(next_nonsingleton->first > cell->first);
	}
      else
	{
	  assert(cell != next_nonsingleton);
	  assert(cell->next_nonsingleton == 0);
	  assert(cell->prev_nonsingleton == 0);
	}
    }
  assert(next_nonsingleton == 0);
  if(do_print)
    fprintf(stderr, "\n");
#endif
}



Partition::BacktrackPoint
Partition::set_backtrack_point()
{
  BacktrackInfo info;
  info.refinement_stack_size = refinement_stack.size();
  BacktrackPoint p = bt_stack.size();
  bt_stack.push_back(info);
  return p;
}



void
Partition::goto_backtrack_point(BacktrackPoint p)
{
  BLISS_ASSERT(p < bt_stack.size());
  BacktrackInfo info = bt_stack[p];
  bt_stack.resize(p);

  const unsigned int dest_refinement_stack_size = info.refinement_stack_size;

  assert(refinement_stack.size() >= dest_refinement_stack_size);
  while(refinement_stack.size() > dest_refinement_stack_size)
    {
      RefInfo i = refinement_stack.pop();
      const unsigned int first = i.split_cell_first;
      Cell *cell = element_to_cell_map[elements[first]];

      if(cell->first != first)
	{
	  assert(cell->split_level <= dest_refinement_stack_size);
	  goto done;
	}
      if(cell->split_level <= dest_refinement_stack_size)
	{
	  goto done;
	}

      {
	const unsigned int new_first = cell->first;
	do
	  {
	    Cell * const next_cell = cell->next;
	    assert(next_cell);
	    /* (Pseudo)free cell */
	    cell->first = 0;
	    cell->length = 0;
	    cell->next->prev_next_ptr = cell->prev_next_ptr;
	    *(cell->prev_next_ptr) = cell->next;
	    cell->next = free_cells;
	    if(cell->next)
	      cell->next->prev_next_ptr = &(cell->next);
	    cell->prev_next_ptr = &free_cells;
	    free_cells = cell;
	    
	    cell = next_cell;
	  }
	while(cell->split_level > dest_refinement_stack_size);

	/* Update element_to_cell_map values of elements added in cell */
	unsigned int *ep = elements + new_first;
	unsigned int * const lp = elements + cell->first;
	while(ep < lp) {
	  element_to_cell_map[*ep] = cell;
	  ep++;
	}
	/* Update cell parameters */
	cell->length = (cell->first + cell->length) - new_first;
	cell->first = new_first;
      }

    done:
      if(i.prev_nonsingleton_first >= 0)
	{
	  Cell * const prev_cell = element_to_cell_map[elements[i.prev_nonsingleton_first]];
	  BLISS_ASSERT(prev_cell->length > 1);
	  cell->prev_nonsingleton = prev_cell;
	  prev_cell->next_nonsingleton = cell;
	}
      else
	{
	  //assert(cell->prev_nonsingleton == 0);
	  cell->prev_nonsingleton = 0;
	  first_nonsingleton_cell = cell;
	}

      if(i.next_nonsingleton_first >= 0)
	{
	  Cell * const next_cell = element_to_cell_map[elements[i.next_nonsingleton_first]];
	  BLISS_ASSERT(next_cell->length > 1);
	  cell->next_nonsingleton = next_cell;
	  next_cell->prev_nonsingleton = cell;
	}
      else
	{
	  //assert(cell->next_nonsingleton == 0);
	  cell->next_nonsingleton = 0;
	}
    }

  consistency_check();
}



Partition::Cell *
Partition::individualize(Partition::Cell * const cell,
			 const unsigned int element)
{
  BLISS_ASSERT(!cell->is_unit());

  unsigned int * const pos = in_pos[element];
  BLISS_ASSERT((unsigned int)(pos - elements) >= cell->first);
  BLISS_ASSERT((unsigned int)(pos - elements) < cell->first + cell->length);
  BLISS_ASSERT(*pos == element);

  const unsigned int last = cell->first + cell->length - 1;
  *pos = elements[last];
  in_pos[*pos] = pos;
  elements[last] = element;
  in_pos[element] = elements + last;
  
  Partition::Cell * const new_cell = aux_split_in_two(cell, cell->length-1);
  BLISS_ASSERT(elements[new_cell->first] == element);
  element_to_cell_map[element] = new_cell;

  consistency_check();

  return new_cell;
} 



Partition::Cell *
Partition::aux_split_in_two(Partition::Cell * const cell,
			    const unsigned int first_half_size)
{
  RefInfo i;

  BLISS_ASSERT(0 < first_half_size && first_half_size < cell->length);

  /* (Pseudo)allocate new cell */
  Cell * const new_cell = free_cells;
  BLISS_ASSERT(new_cell != 0);
  free_cells = new_cell->next;
  /* Update new cell parameters */
  new_cell->first = cell->first + first_half_size;
  new_cell->length = cell->length - first_half_size;
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev_next_ptr = &(new_cell->next);
  new_cell->prev_next_ptr = &(cell->next);
  new_cell->split_level = cell->split_level;
  /* Update old, splitted cell parameters */
  cell->length = first_half_size;
  cell->next = new_cell;
  cell->split_level = refinement_stack.size()+1;
  
  /* Add cell in refinement_stack for backtracking */
  i.split_cell_first = cell->first;
  if(cell->prev_nonsingleton)
    i.prev_nonsingleton_first = cell->prev_nonsingleton->first;
  else
    i.prev_nonsingleton_first = -1;
  if(cell->next_nonsingleton)
    i.next_nonsingleton_first = cell->next_nonsingleton->first;
  else
    i.next_nonsingleton_first = -1;
  refinement_stack.push(i);

  /* Modify nonsingleton cell list */
  if(new_cell->length > 1)
    {
      new_cell->prev_nonsingleton = cell;
      new_cell->next_nonsingleton = cell->next_nonsingleton;
      if(new_cell->next_nonsingleton)
	new_cell->next_nonsingleton->prev_nonsingleton = new_cell;
      cell->next_nonsingleton = new_cell;
    }
  else
    {
      new_cell->next_nonsingleton = 0;
      new_cell->prev_nonsingleton = 0;
    }

  if(cell->is_unit())
    {
      if(cell->prev_nonsingleton)
	cell->prev_nonsingleton->next_nonsingleton = cell->next_nonsingleton;
      else
	first_nonsingleton_cell = cell->next_nonsingleton;
      if(cell->next_nonsingleton)
	cell->next_nonsingleton->prev_nonsingleton = cell->prev_nonsingleton;
      cell->next_nonsingleton = 0;
      cell->prev_nonsingleton = 0;
    }

  return new_cell;
} 



void Partition::print(FILE * const fp)
{
  const char *cell_sep = "";
  fprintf(fp, "[");
  for(Cell *cell = first_cell; cell; cell = cell->next)
    {
      fprintf(fp, "%s", cell_sep);
      cell_sep = ",";
      /* Print cell */
      const char *elem_sep = "";
      fprintf(fp, "{");
      for(unsigned int i = 0; i < cell->length; i++)
	{
	  fprintf(fp, "%s%u", elem_sep, elements[cell->first + i]);
	  elem_sep = ",";
	}
      fprintf(fp, "}");
    }
  fprintf(fp, "]");
}



void Partition::print_signature(FILE * const fp)
{
  const char *cell_sep = "";
  fprintf(fp, "[");
  for(Cell *cell = first_cell; cell; cell = cell->next)
    {
      fprintf(fp, "%s%u", cell_sep, cell->length);
      cell_sep = ",";
    }
  fprintf(fp, "]");
}



void Partition::add_in_splitting_queue(Cell * const cell)
{
  static const unsigned int smallish_cell_threshold = 1;
  assert(!cell->in_splitting_queue);
  cell->in_splitting_queue = true;
  if(cell->length <= smallish_cell_threshold)
    splitting_queue.push_front(cell);
  else
    splitting_queue.push_back(cell);    
}



void Partition::clear_splitting_queue()
{
  while(!splitting_queue.is_empty())
    {
      Cell * const cell = splitting_queue.pop_front();
      assert(cell->in_splitting_queue);
      cell->in_splitting_queue = false;
    }
}



/*
 * Assumes that the invariant values are NOT the same
 * and that the cell contains more than one element
 */
Partition::Cell *
Partition::sort_and_split_cell1(Partition::Cell *cell)
{
#if defined(BLISS_EXPENSIVE_CONSISTENCY_CHECKS)
  assert(cell->length > 1);
  assert(cell->first + cell->length <= N);
  unsigned int nof_0_found = 0;
  unsigned int nof_1_found = 0;
  for(unsigned int i = 0; i < cell->length; i++)
    {
      if(invariant_values[elements[cell->first + i]] == 0)
	nof_0_found++;
      else if(invariant_values[elements[cell->first + i]] == 1)
	nof_1_found++;
      else
	assert(should_not_happen);
    }
  assert(nof_0_found);
  assert(nof_1_found);
  assert(nof_1_found == cell->max_ival_count);
  assert(cell->max_ival == 1);
#endif

  consistency_check();

  /* Allocate new cell */
  Cell *new_cell = free_cells;
  BLISS_ASSERT(new_cell != 0);
  free_cells = new_cell->next;
  if(free_cells) free_cells->prev_next_ptr = &(free_cells);

#define NEW_SORT1
#ifdef NEW_SORT1
  unsigned int *ep0 = elements + cell->first;
  unsigned int *ep1 = ep0 + cell->length - cell->max_ival_count;
  if(cell->max_ival_count > cell->length / 2)
    {
      /* There are more ones than zeros, only move zeros */
      unsigned int * const end = ep0 + cell->length;
      while(ep1 < end)
	{
	  while(invariant_values[*ep1] == 0)
	    {
	      const unsigned int tmp = *ep1;
	      *ep1 = *ep0;
	      *ep0 = tmp;
	      in_pos[tmp] = ep0;
	      in_pos[*ep1] = ep1;
	      ep0++;
	    }
	  element_to_cell_map[*ep1] = new_cell;
	  invariant_values[*ep1] = 0;
	  ep1++;
	}
    }
  else
    {
      /* There are more zeros than ones, only move ones */
      unsigned int * const end = ep1;
      while(ep0 < end)
	{
	  while(invariant_values[*ep0] != 0)
	    {
	      const unsigned int tmp = *ep0;
	      *ep0 = *ep1;
	      *ep1 = tmp;
	      in_pos[tmp] = ep1;
	      in_pos[*ep0] = ep0;
	      ep1++;
	    }
	  ep0++;
	}
      ep1 = end;
      while(ep1 < elements + cell->first + cell->length)
	{
	  element_to_cell_map[*ep1] = new_cell;
	  invariant_values[*ep1] = 0;
	  ep1++;
	}
    }
  /* Update new cell parameters */
  new_cell->first = cell->first + cell->length - cell->max_ival_count;
  new_cell->length = cell->length - (new_cell->first - cell->first);
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev_next_ptr = &(new_cell->next);
  new_cell->prev_next_ptr = &(cell->next);
  new_cell->split_level = cell->split_level;
  /* Update old, splitted cell parameters */
  cell->length = new_cell->first - cell->first;
  cell->next = new_cell;
  cell->split_level = refinement_stack.size()+1;

#else
  /* Sort vertices in the cell according to the invariant values */
  unsigned int *ep0 = elements + cell->first;
  unsigned int *ep1 = ep0 + cell->length;
  while(ep1 > ep0)
    {
      const unsigned int element = *ep0;
      const unsigned int ival = invariant_values[element];
      invariant_values[element] = 0;
      BLISS_ASSERT(ival <= 1);
      BLISS_ASSERT(element_to_cell_map[element] == cell);
      BLISS_ASSERT(in_pos[element] == ep0);
      if(ival == 0)
	{
	  ep0++;
	}
      else
	{
	  ep1--;
	  *ep0 = *ep1;
	  *ep1 = element;
	  element_to_cell_map[element] = new_cell;
	  in_pos[element] = ep1;
	  in_pos[*ep0] = ep0;
	}
    }

  BLISS_ASSERT(ep1 != elements + cell->first);
  BLISS_ASSERT(ep0 != elements + cell->first + cell->length);

  /* Update new cell parameters */
  new_cell->first = ep1 - elements;
  new_cell->length = cell->length - (new_cell->first - cell->first);
  new_cell->next = cell->next;
  if(new_cell->next)
    new_cell->next->prev_next_ptr = &(new_cell->next);
  new_cell->prev_next_ptr = &(cell->next);
  new_cell->split_level = cell->split_level;
  /* Update old, splitted cell parameters */
  cell->length = new_cell->first - cell->first;
  cell->next = new_cell;
  cell->split_level = refinement_stack.size()+1;

#endif /* ifdef NEW_SORT1*/

  /* Add cell in refinement stack for backtracking */
  {
    RefInfo i;
    i.split_cell_first = cell->first;
    if(cell->prev_nonsingleton)
      i.prev_nonsingleton_first = cell->prev_nonsingleton->first;
    else
      i.prev_nonsingleton_first = -1;
    if(cell->next_nonsingleton)
      i.next_nonsingleton_first = cell->next_nonsingleton->first;
    else
      i.next_nonsingleton_first = -1;
    /* Modify nonsingleton cell list */
    if(new_cell->length > 1)
      {
	new_cell->prev_nonsingleton = cell;
	new_cell->next_nonsingleton = cell->next_nonsingleton;
	if(new_cell->next_nonsingleton)
	  new_cell->next_nonsingleton->prev_nonsingleton = new_cell;
	cell->next_nonsingleton = new_cell;
      }
    else
      {
	new_cell->next_nonsingleton = 0;
	new_cell->prev_nonsingleton = 0;
      }
    if(cell->is_unit())
      {
	if(cell->prev_nonsingleton)
	  cell->prev_nonsingleton->next_nonsingleton = cell->next_nonsingleton;
	else
	  first_nonsingleton_cell = cell->next_nonsingleton;
	if(cell->next_nonsingleton)
	  cell->next_nonsingleton->prev_nonsingleton = cell->prev_nonsingleton;
	cell->next_nonsingleton = 0;
	cell->prev_nonsingleton = 0;
      }
    refinement_stack.push(i);
  }


  /* Add cells in splitting queue */
  BLISS_ASSERT(!new_cell->in_splitting_queue);
  if(cell->in_splitting_queue) {
    /* Both cells must be included in splitting_queue in order to have
       refinement to equitable partition */
    add_in_splitting_queue(new_cell);
  } else {
    Cell *min_cell, *max_cell;
    if(cell->length <= new_cell->length) {
      min_cell = cell;
      max_cell = new_cell;
    } else {
      min_cell = new_cell;
      max_cell = cell;
    }
    /* Put the smaller cell in splitting_queue */
    add_in_splitting_queue(min_cell);
    if(max_cell->is_unit()) {
      /* Put the "larger" cell also in splitting_queue */
      add_in_splitting_queue(max_cell);
    }
  }

  consistency_check();

  return new_cell;
}





/**
 * An auxiliary function for distribution count sorting.
 * Build start array so that
 * dcs_start[0] = 0 and dcs_start[i+1] = dcs_start[i] + dcs_count[i].
 */
void Partition::dcs_cumulate_count(const unsigned int max) 
{
  BLISS_ASSERT(max <= 255);
  unsigned int *count_p = dcs_count;
  unsigned int *start_p = dcs_start;
  unsigned int sum = 0;
  for(unsigned int i = max+1; i > 0; i--)
    {
      *start_p = sum;
      start_p++;
      sum += *count_p;
      count_p++;
    }
}


/**
 * Distribution count sorting of cells with invariant values less than 256.
 */
Partition::Cell *
Partition::sort_and_split_cell255(Partition::Cell * const cell,
				  const unsigned int max_ival)
{
  BLISS_ASSERT(max_ival <= 255);

  if(cell->is_unit())
    {
      /* Reset invariant value */
      invariant_values[elements[cell->first]] = 0;
      return cell;
    }
  
#ifdef BLISS_CONSISTENCY_CHECKS
  for(unsigned int i = 0; i < 256; i++)
    assert(dcs_count[i] == 0);
#endif

  /*
   * Compute the distribution of invariant values to the count array
   */
  {
    const unsigned int *ep = elements + cell->first;
    BLISS_ASSERT(element_to_cell_map[*ep] == cell);
    const unsigned int ival = invariant_values[*ep];
    BLISS_ASSERT(ival <= 255);
    dcs_count[ival]++;
    ep++;
#if defined(BLISS_CONSISTENCY_CHECKS)
    bool equal_invariant_values = true;
#endif
    for(unsigned int i = cell->length - 1; i != 0; i--)
      {
	BLISS_ASSERT(element_to_cell_map[*ep] == cell);
	const unsigned int ival2 = invariant_values[*ep];
	BLISS_ASSERT(ival2 <= 255);
	BLISS_ASSERT(ival2 <= max_ival);
	dcs_count[ival2]++;
#if defined(BLISS_CONSISTENCY_CHECKS)
	if(ival2 != ival) {
	  equal_invariant_values = false;
	}
#endif
	ep++;
      }
#if defined(BLISS_CONSISTENCY_CHECKS)
    assert(!equal_invariant_values);
    if(equal_invariant_values) {
      assert(dcs_count[ival] == cell->length);
      dcs_count[ival] = 0;
      clear_ivs(cell);
      return cell;
    }
#endif
  }

  /* Build start array */
  dcs_cumulate_count(max_ival);

  //BLISS_ASSERT(dcs_start[255] + dcs_count[255] == cell->length);
  BLISS_ASSERT(dcs_start[max_ival] + dcs_count[max_ival] == cell->length);

  /* Do the sorting */
  for(unsigned int i = 0; i <= max_ival; i++)
    {
      unsigned int *ep = elements + cell->first + dcs_start[i];
      for(unsigned int j = dcs_count[i]; j > 0; j--)
	{
	  while(true)
	    {
	      const unsigned int element = *ep;
	      const unsigned int ival = invariant_values[element];
	      if(ival == i)
		break;
	      BLISS_ASSERT(ival > i);
	      BLISS_ASSERT(dcs_count[ival] > 0);
	      *ep = elements[cell->first + dcs_start[ival]];
	      elements[cell->first + dcs_start[ival]] = element;
	      dcs_start[ival]++;
	      dcs_count[ival]--;
	    }
	  ep++;
	}
      dcs_count[i] = 0;
    }

#if defined(BLISS_CONSISTENCY_CHECKS)
  for(unsigned int i = 0; i < 256; i++)
    assert(dcs_count[i] == 0);
#endif

  /* split cell */
  Cell * const cell2 = split_cell(cell);
  BLISS_ASSERT(cell2 != cell);
  return cell2;
}


/*
 * Sort the elements in a cell according to their invariant values.
 * The invariant values are not cleared.
 * Warning: the in_pos array is left in incorrect state.
 */
bool Partition::shellsort_cell(Partition::Cell *cell)
{
  unsigned int h;
  unsigned int *ep;

  //BLISS_ASSERT(cell->first + cell->length <= N);

  if(cell->is_unit())
    return false;

  /* Check whether all the elements have the same invariant value */
  bool equal_invariant_values = true;
  {
    ep = elements + cell->first;
    const unsigned int ival = invariant_values[*ep];
    BLISS_ASSERT(element_to_cell_map[*ep] == cell);
    ep++;
    for(unsigned int i = cell->length - 1; i > 0; i--)
      {
	BLISS_ASSERT(element_to_cell_map[*ep] == cell);
	if(invariant_values[*ep] != ival) {
	  equal_invariant_values = false;
	  break;
	}
	ep++;
      }
  }
  if(equal_invariant_values)
    return false;

  ep = elements + cell->first;

  for(h = 1; h <= cell->length/9; h = 3*h + 1)
    ;
  for( ; h > 0; h = h/3) {
    for(unsigned int i = h; i < cell->length; i++) {
      const unsigned int element = ep[i];
      const unsigned int ival = invariant_values[element];
      unsigned int j = i;
      while(j >= h && invariant_values[ep[j-h]] > ival) {
        ep[j] = ep[j-h];
        j -= h;
      }
      ep[j] = element;
    }
  }
  return true;
}


void Partition::clear_ivs(Cell * const cell)
{
  unsigned int *ep = elements + cell->first;
  for(unsigned int i = cell->length; i > 0; i--, ep++)
    invariant_values[*ep] = 0;
}


/*
 * Assumes that the elements in the cell are sorted according to their
 * invariant values.
 */
Partition::Cell *
Partition::split_cell(Partition::Cell * const original_cell)
{
  Cell *cell = original_cell;
  const bool original_cell_was_in_splitting_queue =
    original_cell->in_splitting_queue;
  Cell *largest_new_cell = 0;

  while(true) 
    {
      unsigned int *ep = elements + cell->first;
      const unsigned int * const lp = ep + cell->length;
      const unsigned int ival = invariant_values[*ep];
      invariant_values[*ep] = 0;
      element_to_cell_map[*ep] = cell;
      in_pos[*ep] = ep;
      ep++;
      while(ep < lp)
	{
	  const unsigned int e = *ep;
	  if(invariant_values[e] != ival)
	    break;
	  invariant_values[e] = 0;
	  in_pos[e] = ep;
	  ep++;
	  element_to_cell_map[e] = cell;
	}
      if(ep == lp)
	break;
      
      Cell * const new_cell = aux_split_in_two(cell,
					       (ep - elements) - cell->first);
      
      if(graph->compute_eqref_hash)
	{
	  graph->eqref_hash.update(new_cell->first);
	  graph->eqref_hash.update(new_cell->length);
	  graph->eqref_hash.update(ival);
	}
      
      /* Add cells in splitting_queue */
      assert(!new_cell->in_splitting_queue);
      if(original_cell_was_in_splitting_queue)
	{
	  /* In this case, all new cells are inserted in splitting_queue */
	  assert(cell->in_splitting_queue);
	  add_in_splitting_queue(new_cell);
	}
      else
	{
	  /* Otherwise, we can omit one new cell from splitting_queue */
	  assert(!cell->in_splitting_queue);
	  if(largest_new_cell == 0) {
	    largest_new_cell = cell;
	  } else {
	    assert(!largest_new_cell->in_splitting_queue);
	    if(cell->length > largest_new_cell->length) {
	      add_in_splitting_queue(largest_new_cell);
	      largest_new_cell = cell;
	    } else {
	      add_in_splitting_queue(cell);
	    }
	  }
	}
      /* Process the rest of the cell */
      cell = new_cell;
    }

  consistency_check();
  
  if(original_cell == cell) {
    /* All the elements in cell had the same invariant value */
    return cell;
  }

  /* Add cells in splitting_queue */
  if(!original_cell_was_in_splitting_queue)
    {
      /* Also consider the last new cell */
      assert(largest_new_cell);
      if(cell->length > largest_new_cell->length)
	{
	  add_in_splitting_queue(largest_new_cell);
	  largest_new_cell = cell;
	}
      else
	{
	  add_in_splitting_queue(cell);
	}
      if(largest_new_cell->is_unit())
	{
	  /* Needed in certificate computation */
	  add_in_splitting_queue(largest_new_cell);
	}
    }

  return cell;
}


Partition::Cell *
Partition::zplit_cell(Partition::Cell * const cell,
		      const bool max_ival_info_ok)
{
  BLISS_ASSERT(cell != 0);

  Cell *last_new_cell = cell;

  if(!max_ival_info_ok)
    {
      /* Compute max_ival info */
      assert(cell->max_ival == 0);
      assert(cell->max_ival_count == 0);
      unsigned int *ep = elements + cell->first;
      for(unsigned int i = cell->length; i > 0; i--, ep++)
	{
	  const unsigned int ival = invariant_values[*ep];
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
    }

#ifdef BLISS_CONSISTENCY_CHECKS
  /* Verify max_ival info */
  {
    unsigned int max_ival = 0;
    unsigned int max_ival_count = 0;
    unsigned int *ep = elements + cell->first;
    for(unsigned int i = cell->length; i > 0; i--, ep++)
      {
	const unsigned int ival = invariant_values[*ep];
	if(ival > max_ival)
	  {
	    max_ival = ival;
	    max_ival_count = 1;
	  }
	else if(ival == max_ival)
	  {
	    max_ival_count++;
	  }
      }
    assert(max_ival == cell->max_ival);
    assert(max_ival_count == cell->max_ival_count);
  }
#endif

  /* max_ival info has been computed */

  if(cell->max_ival_count == cell->length)
    {
      /* All invariant values are the same */
      if(cell->max_ival > 0)
	clear_ivs(cell);
      goto done;
    }

  /* All invariant values are not the same */
  if(cell->max_ival == 1)
    {
      /* Specialized splitting for cells with binary invariant values */
      last_new_cell = sort_and_split_cell1(cell);
      goto done;
    }
  if(cell->max_ival < 256)
    {
      /* Specialized splitting for cells with invariant values < 256 */
      last_new_cell = sort_and_split_cell255(cell, cell->max_ival);
      goto done;
    }
  {
    /* Generic sorting and splitting */
    const bool sorted = shellsort_cell(cell);
    assert(sorted);
    last_new_cell = split_cell(cell);
    goto done;
  }

 done:
  cell->max_ival = 0;
  cell->max_ival_count = 0;
  return last_new_cell;
}



} // namespace bliss
