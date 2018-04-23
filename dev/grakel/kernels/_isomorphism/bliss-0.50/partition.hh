#ifndef BLISS_PARTITION_HH
#define BLISS_PARTITION_HH

/*
 * Copyright (c) Tommi Junttila
 * Released under the GNU General Public License version 2.
 */

namespace bliss {
  class Partition;
}

#include <cstdlib>
#include <cstdio>
#include "kstack.hh"
#include "kqueue.hh"
#include "heap.hh"
#include "orbit.hh"
#include "graph.hh"


namespace bliss {

/**
 * \brief A class for refinable, backtrackable ordered partitions.
 *
 * This is rather a data structure with some helper functions than
 * a proper self-contained class.
 * That is, for efficiency reasons the fields of this class are directly
 * manipulated from bliss::AbstractGraph and its subclasses.
 * Conversely, some methods of this class modify the fields of
 * bliss::AbstractGraph, too.
 */
class Partition
{
public:
  /**
   * \brief Data structure for holding information about a cell in a Partition.
   */
  class Cell
  {
  public:
    unsigned int length;
    /* Index of the first element of the cell in
       the Partition::elements array */
    unsigned int first;
    unsigned int max_ival;
    unsigned int max_ival_count;
    bool in_neighbour_heap;
    bool in_splitting_queue;
    /* Pointer to the next cell, null if this is the last one. */
    Cell *next;
    Cell **prev_next_ptr;
    Cell *next_nonsingleton;
    Cell *prev_nonsingleton;
    unsigned int split_level;
    /** Is this a unit cell? */
    bool is_unit() const {return(length == 1); }
  };


private:

  /** \internal
   * Data structure for remembering information about splits in order to
   * perform efficient backtracking over the splits.
   */
  class RefInfo {
  public:
    unsigned int split_cell_first;
    int prev_nonsingleton_first;
    int next_nonsingleton_first;
  };
  /** \internal
   * A stack for remembering the splits, used for backtracking.
   */
  KStack<RefInfo> refinement_stack;

  class BacktrackInfo {
  public:
    unsigned int refinement_stack_size;
  };
  /** \internal
   * The main stack for enabling backtracking.
   */
  std::vector<BacktrackInfo> bt_stack;

public:
  AbstractGraph *graph;

  /* Used during equitable partition refinement */
  KQueue<Cell *> splitting_queue;
  void add_in_splitting_queue(Cell * const cell);
  void clear_splitting_queue();


  /** Type for backtracking points. */
  typedef unsigned int BacktrackPoint;

  /**
   * Get a new backtrack point for the current partition
   */
  BacktrackPoint set_backtrack_point();

  /**
   * Backtrack to the point \a p and remove it.
   */
  void goto_backtrack_point(BacktrackPoint p);

  /**
   * Split the non-unit Cell \a cell = {\a element,e1,e2,...,en} containing
   * the element \a element in two:
   * \a cell = {e1,...,en} and \a newcell = {\a element}.
   * @param cell     a non-unit Cell
   * @param element  an element in \a cell
   * @return         the new unit Cell \a newcell
   */
  Cell *individualize(Cell * const cell,
		      const unsigned int element);

  Cell *aux_split_in_two(Cell * const cell,
			 const unsigned int first_half_size);

  /* The current search level */
  unsigned int level;

  void consistency_check();

private:
  unsigned int N;
  Cell *cells;
  Cell *free_cells;
public:
  Cell *first_cell;
  Cell *first_nonsingleton_cell;
  unsigned int *elements;
  /* invariant_values[e] gives the invariant value of the element e */
  unsigned int *invariant_values;
  /* element_to_cell_map[e] gives the cell of the element e */
  Cell **element_to_cell_map;
  /* in_pos[e] points to the elements array s.t. *in_pos[e] = e  */
  unsigned int **in_pos;

  Partition();
  ~Partition();

  /**
   * Initialize the partition to the unit partition (all elements in one cell)
   * over the \a N > 0 elements {0,...,\a N-1}.
   */
  void init(const unsigned int N);

  /**
   * Returns true iff the partition is discrete, meaning that all
   * the elements are in their own cells.
   */
  bool is_discrete() const {return(free_cells == 0); }

  /**
   * Print the partition into the file stream \a fp.
   */
  void print(FILE * const fp);

  /**
   * Print the partition cell sizes into the file stream \a fp.
   */
  void print_signature(FILE * const fp);

  /*
   * Splits the Cell \a cell into [cell_1,...,cell_n]
   * according to the invariant_values of the elements in \a cell.
   * After splitting, cell_1 == \a cell.
   * Returns the pointer to the Cell cell_n;
   * cell_n != cell iff the Cell \a cell was actually splitted.
   * The flag \a max_ival_info_ok indicates whether the max_ival and
   * max_ival_count fields of the Cell \a cell have consistent values
   * when the method is called.
   * Clears the invariant values of elements in the Cell \a cell as well as
   * the max_ival and max_ival_count fields of the Cell \a cell.
   */
  Cell *zplit_cell(Cell * const cell, const bool max_ival_info_ok);

private:
  /* Auxiliary routines for sorting and splitting cells */

  /** Clear the invariant_values of the elements in the Cell \a cell. */
  void clear_ivs(Cell * const cell);

  Cell *sort_and_split_cell1(Cell *cell);
  Cell *sort_and_split_cell255(Cell * const cell, const unsigned int max_ival);
  bool shellsort_cell(Cell *cell);
  Cell *split_cell(Cell * const cell);

  /*
   * Some auxiliary stuff needed for distribution count sorting.
   * To make the code thread-safe (modulo the requirement that each graph is
   * only accessed in one thread at a time), the arrays are owned by
   * the partition instance, not statically defined.
   */
  unsigned int dcs_count[256];
  unsigned int dcs_start[256];
  void dcs_cumulate_count(const unsigned int max);
};

} // namespace bliss

#endif
