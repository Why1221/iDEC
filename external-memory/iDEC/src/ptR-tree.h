/*****************************************************************
this is the implementation for an R-tree that has an improved storage
scheme for point data. Specifically, now each leaf entry keeps only
d coordinates, as opposed to 2d in an R-tree indexing rectangles.
*****************************************************************/

#ifndef __PTRTREE_H
#define __PTRTREE_H

#include "rtree.h"
#include "rtnode.h"
#include "entry.h"

class ptRTNode;

class ptEntry :public Entry
{
	virtual int get_size();							//gets the number of bytes the entry occupies on disk
	virtual void read_from_buffer(char *_buffer);	//loads the content of an entry from a disk buffer page
	virtual void write_to_buffer(char *_buffer);	//writes the content of an entry to a disk buffer page
	virtual RTNode* new_one_node();				//initiates a node. Override this function to create an object
													//of a class that inherits RTNode
};

class ptRTNode :public RTNode
{
	virtual Entry* new_one_entry();			//initiates an entry. Override this function to create an object
											//of a class that inherits Entry
	virtual RTNode* new_one_node();			//initiates a node. Override this function to create an object
											//of a class that inherits RTNode
};

class ptRTree :public RTree
{
	virtual void fread_next_entry(FILE *_fp, Entry *_d);	//retrieves the next entry from a data file; override this function to
															//load data of different formats
	virtual Entry* new_one_entry();			//initiates an entry. Override this function to create an object
											//of a class that inherits Entry
	virtual RTNode* new_one_node();			//initiates a node. Override this function to create an object
											//of a class that inherits RTNode
};

#endif