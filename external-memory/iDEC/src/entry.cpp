#include <stdio.h>
#include <math.h>
#include <string.h>
#include "entry.h"
#include "rtnode.h"
#include "rtree.h"
#include "linlist.h"
#include "gendef.h"
//------------------------------------------------------------
Entry::Entry()  //this constructor does nothing. call init to do the real work
{
  son_ptr = NULL;
  bounces = NULL;
  my_tree = NULL;
  level = -1;
}
//------------------------------------------------------------
//initialize the basic fields of an entry

void Entry::init(RTree *_rt)
{
  if (_rt)
  {
    dimension = _rt->dimension;
    my_tree = _rt;
    bounces = new float[2 * dimension];
  }
  son = 0;
  //level = 0;	//note that here I deliberately disabled the initialization of "level". this is to facilitate
  //implementing an index where entries of different levels may differ in contents. it is VERY
  //IMPORTANT to set this field by yourself, whenever you initiate an entry. this is necessary
  //for reading/writing the correct information to/from the disk.
  //before reading and writing anything, functions "read_from_buffer" and "write_to_buffer"
  //will first check if level equals -1, in which case the level has not been properly set.
}
//------------------------------------------------------------
Entry::~Entry()
{
  if (bounces)
  {
    delete [] bounces; bounces = NULL;
  }
		
  if (son_ptr != NULL)
  {
    delete son_ptr; son_ptr = NULL;
  }

  my_tree = NULL;
}
//------------------------------------------------------------
void Entry::del_son()
{
  if (son_ptr != NULL)
  {
    delete son_ptr;
    son_ptr = NULL;
  }
}	
//------------------------------------------------------------
//generates a linked-list element from the entry

Linkable* Entry::gen_Linkable()
{
  if (level == -1)
  {
    printf("Entry::gen_Linkable detected that the level of an entry has not been properly set.\n");
    error("The most likely cause of this error is that you forgot to initiate the level of an entry	after creating it. \n", true);
  }
  Linkable *new_link = new Linkable(dimension);
  new_link -> son = son;
  for (int i = 0; i < 2 * dimension; i ++)
    new_link -> bounces[i] = bounces[i];
  new_link -> level = level;

  return new_link;
}
//------------------------------------------------------------
//gets the number of bytes the entry occupies on disk

int Entry::get_size()
{
  return 2 * dimension * sizeof(float) + sizeof(int); //for bounces and son
    
}
//------------------------------------------------------------
//loads the child node of an intermediate entry into memory

RTNode* Entry::get_son()
{
  if (son_ptr == NULL)
  {
    son_ptr = new_one_node();
    son_ptr->init(my_tree, son);
  }

  return son_ptr;
}
//------------------------------------------------------------
//void Entry::init_entry(int _dimension, RTree *_rt)
//{
//	dimension = _dimension;
//    my_tree = _rt;
//    bounces = new float[2 * dimension];
//    son_ptr = NULL;
//    son = 0;
//	level = 0;
//}
//------------------------------------------------------------
//loads the content of an entry from a disk page buffer

void Entry::read_from_buffer(char *_buffer)
{
  if (level == -1)
  {
    printf("Entry::read_from_buffer detected that the level of an entry has not been properly set.\n");
    error("The most likely cause of this error is that you forgot to initiate the level of an entry after creating it. \n", true);
  }

  int i;

  i = 2 * dimension * sizeof(float);
  memcpy(bounces, _buffer, i);

  memcpy(&son, &_buffer[i], sizeof(int));
  i += sizeof(int);
}
//------------------------------------------------------------
//checks the topological relationship between the entry's MBR and the given rectangle

SECTION Entry::section(float *mbr)
{
  bool inside;
  bool overlap;

  overlap = TRUE;
  inside = TRUE;

  for (int i = 0; i < dimension; i++)
  {
    if (mbr[2 * i] > bounces[2 * i + 1] ||  mbr[2 * i + 1] < bounces[2 * i])
      overlap = FALSE;
    if (mbr[2 * i] < bounces[2 * i] ||
        mbr[2 * i + 1] > bounces[2 * i + 1])
      inside = FALSE;
  }
  if (inside)
    return INSIDE;
  else if (overlap)
    return OVERLAP;
  else
    return S_NONE;
}
//------------------------------------------------------------
//sets the content of an entry from a lined-list element

void Entry::set_from_Linkable(Linkable *_link)
{
  son = _link -> son;
  dimension = _link -> dimension;
  memcpy(bounces, _link -> bounces, 2 * dimension * sizeof(float));
  level = _link -> level;

  my_tree = NULL;
  son_ptr = NULL;
}
//------------------------------------------------------------
// writes the content of an entry to a disk buffer page

void Entry::write_to_buffer(char *_buffer)
{
  if (level == -1)
  {
    printf("Entry::write_to_buffer detected that the level of an entry has not been properly set.\n");
    error("The most likely cause of this error is that you forgot to initiate the level of an entry after creating it. \n", true);
  }

  int i;

  i = 2 * dimension * sizeof(float);
  memcpy(_buffer, bounces, i);
  memcpy(&_buffer[i], &son, sizeof(int));
  i += sizeof(int);
}
//------------------------------------------------------------
//this function compares two entries based on son, dimension, and extents

bool Entry::check_equal(Entry *_d)
{
  if (son != _d->son) return false;
  if (dimension != _d->dimension) return false;
  for (int i = 0; i < 2 * dimension; i++)
    if (fabs(bounces[i] - _d->bounces[i]) > FLOATZERO) return false;
  return true;
}
//------------------------------------------------------------
//this function assigns all fields of an entry to those of the given entry

void Entry::set_equal_to(Entry *_d)
{
  dimension = _d->dimension;
  son = _d->son;
  son_ptr = _d->son_ptr;
  memcpy(bounces, _d->bounces, sizeof(float) * 2 * dimension);
  my_tree = _d->my_tree;
  level = _d->level;
}

/*****************************************************************
initiates a node. Override this function to create an object of a class that inherits RTNode

Coded by Yufei Tao 05/31/06
*****************************************************************/

RTNode* Entry::new_one_node()
{
  RTNode *nd = new RTNode();
  return nd;
}


