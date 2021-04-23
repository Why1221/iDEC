#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "ptR-tree.h"
#include "gendef.h"
#include <string.h>

//get the number of bytes occupied by an entry on disk.
//coded by Yufei Tao, may 31, 2006.

int ptEntry::get_size()
{
	if (level == -1)
	{
		printf("ptEntry::get_size detected that the level of an entry has not been properly set.\n");
		error("The most likely cause of this error is that you forgot to initiate the level of an entry after creating it. \n", true);
	}

	if (level == 0)
	{
		return dimension * sizeof(float) + sizeof(int); //for bounces and son
	}

	return Entry::get_size();
}

//loads the content of an entry from a disk page buffer.
//coded by Yufei Tao, may 31, 2006.

void ptEntry::read_from_buffer(char *_buffer)
{
	if (level == -1)
	{
		printf("ptEntry::read_from_buffer detected that the level of an entry has not been properly set.\n");
		error("The most likely cause of this error is that you forgot to initiate the level of an entry after creating it. \n", true);
	}

	if (level == 0)
	{
		int i = 0;

		for (int j = 0; j < dimension; j ++)
		{
			memcpy(&bounces[2 * j], &_buffer[i], sizeof(float));
			i += sizeof(float);
			bounces[2 * j + 1] = bounces[2 * j];
		}

		memcpy(&son, &_buffer[i], sizeof(int));
		i += sizeof(int);
	}
	else
		Entry::read_from_buffer(_buffer);
}

// writes the content of an entry to a disk buffer page
//coded by Yufei Tao, may 31, 2006.

void ptEntry::write_to_buffer(char *_buffer)
{
	if (level == -1)
	{
		printf("ptEntry::write_to_buffer detected that the level of an entry has not been properly set.\n");
		error("The most likely cause of this error is that you forgot to initiate the level of an entry after creating it. \n", true);
	}

    if (level == 0)
	{
		int i = 0;

		for (int j = 0; j < dimension; j ++)
		{
			memcpy(&_buffer[i], &bounces[2 * j], sizeof(float));
			i += sizeof(float);
		}

		memcpy(&_buffer[i], &son, sizeof(int));
		i += sizeof(int);
	}
	else
		Entry::write_to_buffer(_buffer);
}

//initiates a node
//Coded by Yufei Tao 05/31/06

RTNode* ptEntry::new_one_node()
{
	ptRTNode *nd = new ptRTNode();
	return nd;
}

//initiates an entry. Override this function to create an object of a class that inherits Entry
//Coded by Yufei Tao 05/31/06

Entry* ptRTNode::new_one_entry()
{
	ptEntry *e = new ptEntry();
	return e;
}

//initiates a node. Override this function to create an object of a class that inherits RTNode
//Coded by Yufei Tao 05/31/06

RTNode* ptRTNode::new_one_node()
{
	ptRTNode *nd = new ptRTNode();
	return nd;
}

//initiates a node. Override this function to create an object of a class that inherits RTNode
//Coded by Yufei Tao 05/31/06

RTNode* ptRTree::new_one_node()
{
	ptRTNode *nd = new ptRTNode();
	return nd;
}

//initiates an entry. Override this function to create an object of a class that inherits Entry
//Coded by Yufei Tao 05/31/06

Entry* ptRTree::new_one_entry()
{
	ptEntry *e = new ptEntry();
	return e;
}

//retrieves the next entry from a data file
//Coded by Yufei Tao 05/31/06

void ptRTree::fread_next_entry(FILE *_fp, Entry *_d)
{
	fscanf(_fp, "%d\t", &(_d->son));
	for (int i = 0; i < dimension; i++)
	{
    	fscanf(_fp, "%f", &(_d->bounces[2*i]));
		_d->bounces[2*i+1] = _d->bounces[2*i];
	}
	fscanf(_fp, "\n");
}
