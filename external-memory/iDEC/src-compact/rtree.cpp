/*rtree.cpp
  this file implements the RTree class*/
#include "rtree.h"
#include <math.h>
#include <string.h>
#include <queue>
#include "blk_file.h"
#include "cache.h"
#include "entry.h"
#include "gendef.h"
#include "heap.h"
#include "linlist.h"
#include "rtnode.h"
//------------------------------------------------------------
RTree::RTree() {
  root_ptr = NULL;
  re_level = NULL;
  re_data_cands = NULL;
  deletelist = NULL;
}
//------------------------------------------------------------
// initializes a new tree

void RTree::init(char *_fname, int _b_length, Cache *_c, int _dimension) {
  FILE *fp = fopen(_fname, "r");
  if (fp) {
    fclose(fp);
    printf("The file \"%s\" exists. Replace? (y/n)", _fname);
    char c = getchar();
    getchar();
    if (c != 'y' && c != 'Y') error("", true);
    remove(_fname);
  }

  file = new BlockFile(_fname, _b_length);
  cache = _c;

  re_data_cands = new LinList();
  deletelist = new LinList();

  dimension = _dimension;
  root = 0;
  root_ptr = NULL;
  root_is_data = TRUE;
  num_of_data = num_of_inodes = num_of_dnodes = 0;

  root_ptr = new_one_node();
  root_ptr->init(0, this);
  num_of_dnodes++;
  root_ptr->level = 0;
  root = root_ptr->block;
  del_root();
}
//------------------------------------------------------------
// loads an existing tree

void RTree::init(char *_fname, Cache *_c) {
  //    int j;

  file = new BlockFile(_fname, 0);
  cache = _c;

  re_data_cands = new LinList();
  deletelist = new LinList();

  char *header = new char[file->get_blocklength()];
  file->read_header(header);
  read_header(header);
  delete[] header;

  root_ptr = NULL;
}
//------------------------------------------------------------
// construct a tree from a set of rectangles

void RTree::build_from_file(char *_dsname) {
  Entry *d;
  FILE *fp;

  int record_count = 0;

  if ((fp = fopen(_dsname, "r")) == NULL) {
    delete this;
    error("Could not open the data file", TRUE);
  } else {
    // int id = 0;
    //	  float x0, y0, x1, y1;
    while (!feof(fp)) {
      record_count++;
      // id ++; //disable this variable if the id comes with the dataset

      // debugging code..............
      // if (record_count == 2011)
      // printf("testing...");
      // printf("%d\n", record_count);
      if (record_count % 1000 == 0) {
        printf("inserting object %d\n", record_count);
      }

      d = new_one_entry();
      d->init(this);
      d->level = 0;

      fread_next_entry(fp, d);

      insert(d);  // d will be deleted in insert()
    }
  }

  fclose(fp);
  printf("\n");

  del_root();
}
//------------------------------------------------------------
void RTree::close() {
  char *header = new char[file->get_blocklength()];
  write_header(header);
  file->set_header(header);
  delete[] header;

  if (root_ptr != NULL) {
    delete root_ptr;
    root_ptr = NULL;
  }

  if (cache) {
    cache->flush();  // note that a cache object is not destroyed here. Remember
                     // that this object
    // was not initiated in the R-tree, but passed in when the tree was
    // initiated.
  }

  delete file;
  file = NULL;
  delete re_data_cands;
  re_data_cands = NULL;
  delete deletelist;
  deletelist = NULL;
}
//------------------------------------------------------------
RTree::~RTree() {
  if (file) {
    printf("The pointer 'file' is not NULL\n");
    error(
        "The most likely cause of this error is that you forgot to call "
        ".close before deleting the object",
        true);
  }

  //	if (cache)
  //	{
  //		printf("The pointer 'cache' is not NULL\n");
  //		error("The most likely cause of this error is that you forgot to
  // call .close before deleting the object", 			true);
  //	}

  if (root_ptr) {
    printf("The pointer 'root_ptr' is not NULL\n");
    error(
        "The most likely cause of this error is that you forgot to call "
        ".close before deleting the object",
        true);
  }

  if (re_data_cands) {
    printf("The pointer 're_data_cands' is not NULL\n");
    error(
        "The most likely cause of this error is that you forgot to call "
        ".close before deleting the object",
        true);
  }

  if (deletelist) {
    printf("The pointer 'deletelist' is not NULL\n");
    error(
        "The most likely cause of this error is that you forgot to call "
        ".close before deleting the object",
        true);
  }

  // printf("This R-Tree contains %d internal, %d data node(s) and %d
  // object(s)\n", num_of_inodes, num_of_dnodes,
  //        num_of_data);
}
//------------------------------------------------------------
// deletes an object

bool RTree::delete_entry(Entry *_d) {
  load_root();

  R_DELETE del_ret;
  del_ret = root_ptr->delete_entry(_d);

  if (del_ret == NOTFOUND) return false;
  if (del_ret == ERASED)
    error("RTree::delete_entry: The root has been deleted\n", true);

  if (root_ptr->level > 0 && root_ptr->num_entries == 1)
  // there is only one entry in the root but the root
  // is not leaf.  in this case, the child of the root is exhalted to root
  {
    root = root_ptr->entries[0]->son;
    delete root_ptr;
    root_ptr = NULL;
    load_root();
    num_of_inodes--;
  }

  num_of_data--;

  // will reinsert some entries
  while (deletelist->get_num() > 0) {
    Linkable *e;
    e = deletelist->get_first();
    Entry *new_e = new_one_entry();
    new_e->init(this);
    new_e->set_from_Linkable(e);
    deletelist->erase();
    insert(new_e);
    num_of_data--;
  }

  delete root_ptr;
  root_ptr = NULL;

  return true;
}
//------------------------------------------------------------
// sees if an object exists

bool RTree::FindLeaf(Entry *_e) {
  load_root();
  bool ret = root_ptr->FindLeaf(_e);
  del_root();
  return ret;
}
//------------------------------------------------------------
// inserts an object

void RTree::insert(Entry *_d) {
  int i, j;
  RTNode *sn;
  RTNode *nroot_ptr;
  int nroot;
  Entry *de;
  R_OVERFLOW split_root;
  Entry *dc;
  // float *nmbr;
  int16_t *nmbr;

  // load root into memory
  load_root();

  // no overflow occured until now
  re_level = new bool[root_ptr->level + 1];
  for (i = 0; i <= root_ptr->level; i++) re_level[i] = FALSE;

  // insert d into re_data_cands as the first entry to insert
  // make a copy of d because it should be erased later
  Linkable *new_link;
  new_link = _d->gen_Linkable();
  re_data_cands->insert(new_link);

  delete _d;  // we follow the convention that the entry will be deleted when
              // insertion finishes

  j = -1;
  while (re_data_cands->get_num() > 0) {
    // first try to insert data, then directory entries
    Linkable *d_cand;
    d_cand = re_data_cands->get_first();
    if (d_cand != NULL) {
      // since "erase" deletes the data itself from the
      // list, we should make a copy of the data before
      // erasing it
      dc = new_one_entry();
      dc->init(this);
      dc->set_from_Linkable(d_cand);
      re_data_cands->erase();

      // start recursive insert with root
      split_root = root_ptr->insert(dc, &sn);
    } else
      error("RTree::insert: inconsistent list re_data_cands", TRUE);

    if (split_root == SPLIT)
    // insert has lead to split --> new root-page with two sons (i.e. root and
    // sn)
    {
      nroot_ptr = new_one_node();
      nroot_ptr->init(root_ptr->level + 1, this);
      num_of_inodes++;
      nroot = nroot_ptr->block;

      de = new_one_entry();
      de->init(this);
      de->level = nroot_ptr->level;
      nmbr = root_ptr->get_mbr();
      // memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
      memcpy(de->bounces, nmbr, 2 * dimension * sizeof(int16_t));
      delete[] nmbr;
      de->son = root_ptr->block;
      de->son_ptr = root_ptr;
      nroot_ptr->enter(de);

      de = new_one_entry();
      de->init(this);
      de->level = nroot_ptr->level;
      nmbr = sn->get_mbr();
      // memcpy(de -> bounces, nmbr, 2*dimension*sizeof(float));
      memcpy(de->bounces, nmbr, 2 * dimension * sizeof(int16_t));
      delete[] nmbr;
      de->son = sn->block;
      de->son_ptr = sn;
      nroot_ptr->enter(de);

      root = nroot;
      root_ptr = nroot_ptr;

      root_is_data = FALSE;
    }
    j++;
  }

  num_of_data++;

  delete[] re_level;

  delete root_ptr;
  root_ptr = NULL;
}
//------------------------------------------------------------
// loads the root into the memory

void RTree::load_root() {
  if (root_ptr == NULL) {
    root_ptr = new_one_node();
    root_ptr->init(this, root);
  }
}
//------------------------------------------------------------
// void RTree::NNQuery(float *QueryPoint,
//					SortedLinList *res)
//{
//      float nearest_distanz;
//
//      // load root node into main memory
//      load_root();
//
//      nearest_distanz = MAXREAL;
//
//      root_ptr->NNSearch(QueryPoint,res,&nearest_distanz);
//
//	  delete root_ptr;
//	  root_ptr = NULL;
//}
//------------------------------------------------------------
// performs a range query

// int RTree::rangeQuery(float *_mbr, SortedLinList *_res, bool _ids_wanted)
int RTree::rangeQuery(int16_t *_mbr, SortedLinList *_res, bool _ids_wanted) {
  load_root();
  int ret = root_ptr->rangeQuery(_mbr, _res, _ids_wanted);
  del_root();
  return ret;
}
//------------------------------------------------------------
// loads the parameters of the tree from the disk buffer page

void RTree::read_header(char *_buffer) {
  int i;

  memcpy(&dimension, _buffer, sizeof(dimension));
  i = sizeof(dimension);

  memcpy(&num_of_data, &_buffer[i], sizeof(num_of_data));
  i += sizeof(num_of_data);

  memcpy(&num_of_dnodes, &_buffer[i], sizeof(num_of_dnodes));
  i += sizeof(num_of_dnodes);

  memcpy(&num_of_inodes, &_buffer[i], sizeof(num_of_inodes));
  i += sizeof(num_of_inodes);

  memcpy(&root_is_data, &_buffer[i], sizeof(root_is_data));
  i += sizeof(root_is_data);

  memcpy(&root, &_buffer[i], sizeof(root));
  i += sizeof(root);
}
//------------------------------------------------------------
// writes the parameters of the tree into the disk buffer page

void RTree::write_header(char *buffer) {
  int i;

  memcpy(buffer, &dimension, sizeof(dimension));
  i = sizeof(dimension);

  memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
  i += sizeof(num_of_data);

  memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
  i += sizeof(num_of_dnodes);

  memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
  i += sizeof(num_of_inodes);

  memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
  i += sizeof(root_is_data);

  memcpy(&buffer[i], &root, sizeof(root));
  i += sizeof(root);
}

/*****************************************************************
performs a kNN query using the best-first algorithm
para:
q: the query point
k: the numbers of NNs retrieved
hp: the heap for the best-first search
fardist: the distance from q to the kth NN

Coded by Yufei Tao 16/12/02
*****************************************************************/
/*
// int RTree::kNN(float *_q, int _k, Heap *_hp, float &_fardist)
int RTree::kNN(int16_t *_q, int _k, Heap *_hp, float &_fardist) {
  _fardist = MAXREAL;
  int found_cnt = 0;
  int io_cnt = 0;

  // clear the node access summary-----------------------------------
  for (int i = 0; i < 10; i++) na[i] = 0;
  //----------------------------------------------------------------

  // first insert the root into the heap-------------------------
  HeapEntry *he = new HeapEntry();
  he->son1 = root;
  he->key = 0;
  he->level = 1;  // this is not important but make sure it is greather than 0
  _hp->insert(he);
  delete he;
  //------------------------------------------------------------

  while (_hp->used > 0) {
    // remove the top heap-------------------------------------
    HeapEntry *he = new HeapEntry();
    _hp->remove(he);
    int son = he->son1;
    int level = he->level;
    if (level == 0) _fardist = he->key;
    delete he;
    //--------------------------------------------------------

    if (level == 0) {
      //			_rslt[found_cnt]=son;
      found_cnt++;
      if (found_cnt == _k) return io_cnt;
    } else {
      io_cnt++;
      RTNode *child = new_one_node();
      child->init(this, son);
      for (int i = 0; i < child->num_entries; i++) {
        // now init a new heap entry-----------------------
        HeapEntry *he = new HeapEntry();
        he->son1 = child->entries[i]->son;
        he->level = child->level;
        he->key = MINDIST(_q, child->entries[i]->bounces, dimension);
        _hp->insert(he);
        delete he;
        //------------------------------------------------
      }
    
      na[(int)child->level]++;
      delete child;
    }
  }
  // this line should not be reached
  perror("Something goes wrong");
  return io_cnt;
}
*/
// int RTree::kNN(float *_q, int _k, Heap *_hp, float &_fardist)
int RTree::kNN(int16_t *_q, int _k, std::vector<std::pair<size_t, float > >& results) {
  //_fardist = MAXREAL;
  int found_cnt = 0;
  int io_cnt = 0;

  // clear the node access summary-----------------------------------
  for (int i = 0; i < 10; i++) na[i] = 0;
  //----------------------------------------------------------------

  // Using lambda to compare elements.
  auto cmp = [](HeapEntry left, HeapEntry right) { return left.key > right.key; };
  std::priority_queue<HeapEntry, std::vector<HeapEntry>, decltype(cmp)> sq(cmp);
  // first insert the root into the heap-------------------------
  HeapEntry *he = new HeapEntry();
  he->son1 = root;
  he->key = 0;
  he->level = 1;  // this is not important but make sure it is greather than 0
  //_hp->insert(he);
  sq.push(*he);
  delete he;

  //------------------------------------------------------------

  while (!sq.empty()) {
    // remove the top heap-------------------------------------
    //HeapEntry *he = new HeapEntry();
    auto he = sq.top();
    sq.pop();
    int son = he.son1;
    int level = he.level;
    float dist = he.key;
    //--------------------------------------------------------

    if (level == 0) {
      //			_rslt[found_cnt]=son;
      found_cnt++;
      results.emplace_back(son, dist);
      if (found_cnt == _k) return io_cnt;
    } else {
      io_cnt++;
      RTNode *child = new_one_node();
      child->init(this, son);
      for (int i = 0; i < child->num_entries; i++) {
        // now init a new heap entry-----------------------
        HeapEntry *he = new HeapEntry();
        he->son1 = child->entries[i]->son;
        he->level = child->level;
        he->key = MINDIST(_q, child->entries[i]->bounces, dimension);
        sq.push(*he);
        delete he;
        //------------------------------------------------
      }

      na[(int)child->level]++;
      delete child;
    }
  }
  // this line should not be reached
  perror("Something goes wrong");
  return io_cnt;
}
//------------------------------------------------------------
// frees the root from memory

void RTree::del_root() {
  delete root_ptr;
  root_ptr = NULL;
}

/*****************************************************************
initiates a node. Override this function to create an object of a class that
inherits RTNode

Coded by Yufei Tao 05/31/06
*****************************************************************/

RTNode *RTree::new_one_node() {
  RTNode *nd = new RTNode();
  return nd;
}

/*****************************************************************
initiates an entry. Override this function to create an object of a class that
inherits Entry

Coded by Yufei Tao 05/31/06
*****************************************************************/

Entry *RTree::new_one_entry() {
  Entry *e = new Entry();
  return e;
}

/*****************************************************************
retrieves the next entry from a data file

Coded by Yufei Tao 05/31/06
*****************************************************************/

void RTree::fread_next_entry(FILE *_fp, Entry *_d) {
  if (fscanf(_fp, "%d ", &(_d->son)) != 1) {
    perror("fread_next_entry failed");
  }
  for (int i = 0; i < dimension; i++) {
    // fscanf(_fp, "%f %f", &(_d->bounces[2 * i]), &(_d->bounces[2 * i + 1]));
    if (fscanf(_fp, "%hd %hd", &(_d->bounces[2 * i]),&(_d->bounces[2 * i + 1])) != 2) {
      perror("fread_next_entry failed");
    }
#ifdef DEBUG
    printf("%f %f ", _d->bounces[2 * i], _d->bounces[2 * i + 1]);
#endif
  }
  if (fscanf(_fp, "\n") != 0) {
    perror("file format is incorrect");
  }
#ifdef DEBUG
  printf("\n");
#endif
}
