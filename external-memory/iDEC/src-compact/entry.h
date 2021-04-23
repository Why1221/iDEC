#ifndef __Entry
#define __Entry
//------------------------------------------------------------
#include "rtree.h"
#include <cstdint>
//------------------------------------------------------------
class RTNode;
class RTree;
struct Linkable;
//------------------------------------------------------------
class Entry {
public:
  //--===on disk===--
  int son;
  // float *bounces;
  int16_t *bounces;
  //--===others===--
  int dimension;
  int level;
  RTree *my_tree;
  RTNode *son_ptr;

  //--===functions===--
  Entry();
  RTNode *get_son(); // loads the child node of an intermediate entry into
                     // memory
  SECTION section(int16_t *mbr); // checks the topological relationship between
                               // the entry's MBR and the given rectangle

  virtual ~Entry();
  virtual bool check_equal(Entry *_d); // this function compares two entries
                                       // based on son, dimension, and extents
  virtual void
  del_son(); // frees the subtree of an intermediate entry from memory
  virtual int get_size(); // gets the number of bytes the entry occupies on disk
  virtual void init(RTree *_rt); // initialize the basic fields of an entry
  virtual Linkable *
  gen_Linkable(); // generates a linked-list element from the entry. the element
                  // may need to be redefined in implementing an alternative
                  // index.
  virtual RTNode *
  new_one_node(); // initiates a node. Override this function to create an
                  // object of a class that inherits RTNode
  virtual void read_from_buffer(
      char *_buffer); // loads the content of an entry from a disk buffer page
  virtual void set_equal_to(Entry *_d); // this function assigns all fields of
                                        // an entry to those of the given entry
  virtual void set_from_Linkable(
      Linkable
          *_link); // sets the content of an entry from a lined-list element
  virtual void write_to_buffer(
      char *_buffer); // writes the content of an entry to a disk buffer page

  // bool section_circle(float *center, float radius);	//checks
  // void init_entry(int _dimension, RTree *_rt);
};

#endif
