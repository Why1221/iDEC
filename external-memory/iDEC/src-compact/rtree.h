/* rtree.h
   this file defines the class RTree*/

#ifndef __RTREE
#define __RTREE
//------------------------------------------------------------
#include <cstdint>
#include <utility>
#include <vector>
#include "cache.h"
#include "heap.h"
//------------------------------------------------------------
class LinList;
class SortedLinList;
class Cache;
class RTNode;
class Entry;
//------------------------------------------------------------
enum SECTION { OVERLAP, INSIDE, S_NONE };
enum R_OVERFLOW { SPLIT, REINSERT, NONE };
enum R_DELETE { NOTFOUND, NORMAL, ERASED };

class RTree : public Cacheable {
public:
  //--===on disk===--
  int dimension;
  int num_of_data;
  int num_of_dnodes;
  int num_of_inodes;
  int root;
  bool root_is_data;
  //--===others===--
  RTNode *root_ptr;
  bool *re_level;
  LinList *re_data_cands;
  LinList *deletelist;
  //--===added for cost summary===--
  int na[10];

  RTree();

  virtual ~RTree();
  virtual void
  build_from_file(char *_dsname); // construct a tree from a set of rectangles
  virtual void close(); // call this function before destroying the object
  virtual bool delete_entry(Entry *_d); // delets an object
  virtual void del_root();              // frees the root from memory
  virtual bool FindLeaf(Entry *_e);     // sees if an object exists
  virtual void fread_next_entry(
      FILE *_fp,
      Entry *_d); // retrieves the next entry from a data file; override this
                  // function to load data of different formats
  virtual void init(char *_fname, int _b_length, Cache *_c,
                    int _dimension);          // creates a new tree
  virtual void init(char *_fname, Cache *_c); // loads an existing tree
  virtual void insert(Entry *_d);
  /*
  virtual int kNN(int16_t *_q, int _k, Heap *_hp,
                  float &_fardist); // performs a kNN query
  */
  virtual int kNN(int16_t *_q, int _k, std::vector<std::pair<size_t, float > >& results); // performs a kNN query
  virtual void load_root();         // loads the root into memory
  virtual Entry *
  new_one_entry(); // initiates an entry. Override this function to create an
                   // object of a class that inherits Entry
  virtual RTNode *
  new_one_node(); // initiates a node. Override this function to create an
                  // object of a class that inherits RTNode
  virtual int rangeQuery(int16_t *_mbr, SortedLinList *_res,
                         bool _ids_wanted); // performs a range query
  virtual void read_header(char *_buffer);  // loads the parameters of the tree
                                           // from the disk buffer page
  virtual void write_header(char *_buffer); // writes the parameters of the tree
                                            // into the disk buffer page

  // int get_num() { return num_of_data; }
  // void NNQuery(float *QueryPoint, SortedLinList *res);
};

#endif // __RTREE
