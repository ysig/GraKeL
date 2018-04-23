#ifndef BLISS_HEAP_HH
#define BLISS_HEAP_HH

namespace bliss {

/**
 * \brief A capacity bounded heap data structure.
 */

class Heap
{
  unsigned int N;
  unsigned int n;
  unsigned int *array;
  void upheap(unsigned int k);
  void downheap(unsigned int k);
public:
  /**
   * Create a new heap.
   * init() must be called after this.
   */
  Heap() {array = 0; n = 0; N = 0; }
  ~Heap();

  /**
   * Initialize the heap to have the capacity to hold \e size elements.
   */
  void init(const unsigned int size);

  /**
   * Is the heap empty?
   * Time complexity is O(1).
   */
  bool is_empty() const {return(n==0); }

  /**
   * Remove all the elements in the heap.
   * Time complexity is O(1).
   */
  void clear() {n = 0;}

  /**
   * Insert the element \a e in the heap.
   * Time complexity is O(log(N)), where N is the number of elements
   * currently in the heap.
   */
  void insert(const unsigned int e);

  /**
   * Remove and return the smallest element in the heap.
   * Time complexity is O(log(N)), where N is the number of elements
   * currently in the heap.
   */
  unsigned int remove();
};

} // namespace bliss

#endif
